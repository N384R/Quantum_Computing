'''
Unitary Coupled Cluster Singles and Doubles (UCCSD) ansatz group

The class has the following methods:

generate_coeff: Generates UCCSD coefficients.
ansatz: Generates UCCSD ansatz circuit.
'''
from itertools import product, combinations
import numpy as np
import scipy.optimize as opt
from jqc.mapper.jordan_wigner import JordanWignerMapper

def sign_p(val):
    'Sign function'
    return '+' if val > 0 else '-'

def sign_m(val):
    'Sign function'
    return '-' if val > 0 else '+'

def singles(val, i, j):
    'Singles operator'
    operator = f'{sign_p(val)} {abs(val):f} {j}^ {i} ' + '\n'
    operator += f'{sign_m(val)} {abs(val):f} {i}^ {j} ' + '\n'
    return operator

def doubles(val, i, j, k, l):
    'Doubles operator'
    operator = f'{sign_p(val)} {abs(val):f} {k}^ {l}^ {j} {i} ' + '\n'
    operator += f'{sign_m(val)} {abs(val):f} {i}^ {j}^ {l} {k} ' + '\n'
    return operator

def boundary(coeff):
    'Boundary condition'
    return [(-2*np.pi, 2*np.pi)] * len(coeff)

class UCCSD:
    'Unitary Coupled Cluster Singles and Doubles (UCCSD) ansatz'

    @staticmethod
    def call_optimizer(func, coeff, method) -> opt.OptimizeResult:
        'Optimize the coefficients'
        return opt.minimize(func, coeff, method=method, bounds=boundary(coeff))

    def generate_coeff(self, profile, coeff=0.0):
        'Generate UCCSD coefficients'
        no = profile.num_orb
        ne = profile.num_elec // 2
        count = 0
        for _ in range(ne):
            for _ in range(ne, no):
                count += 2

        for _ in product(range(ne), repeat = 2):
            for _ in product(range(ne, no), repeat = 2):
                count += 1

        return [coeff] * count

    def mapping(self, profile, coeff):
        'Generate UCCSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        ne = profile.num_elec // 2
        value = iter(coeff)
        for i in range(ne):
            for j in range(ne, no):
                uccsd_fermion += singles(next(value), i, j)
                uccsd_fermion += singles(next(value), i+no, j+no)

        for i, j in product(range(ne), repeat = 2):
            for k, l in product(range(ne, no), repeat = 2):
                uccsd_fermion += doubles(next(value), i, j+no, k, l+no)

        return JordanWignerMapper(uccsd_fermion)

    def ansatz(self, qc, profile, coeff):
        'Generate UCCSD ansatz circuit'
        for p_string, values in self.mapping(profile, coeff).items():
            chk = []
            for idx, p in p_string.items():
                if p.symbol == 'X':
                    qc.h(idx)
                elif p.symbol == 'Y':
                    qc.sdg(idx)
                    qc.h(idx)
                if p.symbol != 'I':
                    chk.append(idx)

            for i, j in zip(chk, chk[1:]):
                qc.cx(i, j)

            qc.rz(values.imag, max(chk))

            for i, j in zip(chk[:-1][::-1], chk[1:][::-1]):
                qc.cx(i, j)

            for idx, p in p_string.items():
                if p.symbol == 'X':
                    qc.h(idx)
                elif p.symbol == 'Y':
                    qc.h(idx)
                    qc.s(idx)

class eUCCSD(UCCSD):
    ' entangled Unitary Coupled Cluster Singles and Doubles (cUCCSD) ansatz'

    def generate_coeff(self, profile, coeff=0.0):
        'Generate Spin Flip UCCSD coefficients'
        no = profile.num_orb
        ne = profile.num_elec // 2
        count = 0
        for _ in range(ne):
            for _ in range(ne, no):
                count += 2

        for _ in product(range(ne), repeat = 2):
            for _ in product(range(ne, no), repeat = 2):
                count += 1

        for _ in combinations(range(no), 2):
            count += 1

        return [coeff] * count

    def mapping(self, profile, coeff):
        'Generate Spin Flip UCCSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        ne = profile.num_elec // 2
        value = iter(coeff)
        for i in range(ne):
            for j in range(ne, no):
                uccsd_fermion += singles(next(value), i, j)
                uccsd_fermion += singles(next(value), i+no, j+no)

        for i, j in product(range(ne), repeat = 2):
            for k, l in product(range(ne, no), repeat = 2):
                uccsd_fermion += doubles(next(value), i, j+no, k, l+no)

        for i, j in combinations(range(no), 2):
            uccsd_fermion += doubles(next(value), j, i+no, j+no, i)

        return JordanWignerMapper(uccsd_fermion)

class UCCGSD(UCCSD):
    'Unitary Coupled Cluster Generalized Singles and Doubles (UCCGSD) ansatz'

    def generate_coeff(self, profile, coeff=0.0):
        'Generate UCCGSD coefficients'
        no = profile.num_orb
        count = 0
        for _ in combinations(range(no), 2):
            count += 2

        for _ in combinations(range(no), 4):
            count += 2

        for i, j in product(range(no), repeat = 2):
            for _ in product(range(i+1, no), range(j+1, no)):
                count += 1

        return [coeff] * count

    def mapping(self, profile, coeff):
        'Generate UCCGSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        value = iter(coeff)
        for i, j in combinations(range(no), 2):
            uccsd_fermion += singles(next(value), i, j)

        for i, j, k, l in combinations(range(no), 4):
            uccsd_fermion += doubles(next(value), i, j, k, l)
            uccsd_fermion += doubles(next(value), i+no, j+no, k+no, l+no)

        for i, j in product(range(no), repeat = 2):
            for k, l in product(range(i+1, no), range(j+1, no)):
                uccsd_fermion += doubles(next(value), i, j+no, k, l+no)

        return JordanWignerMapper(uccsd_fermion)

class eUCCGSD(UCCSD):
    'entangled Unitary Coupled Cluster Generalized Singles and Doubles (eUCCGSD) ansatz'

    def generate_coeff(self, profile, coeff=0.0):
        'Generate UCCGSD coefficients'
        no = profile.num_orb
        count = 0
        for _ in combinations(range(no), 2):
            count += 2

        for _ in combinations(range(no), 4):
            count += 2

        for i, j in product(range(no), repeat = 2):
            for _ in product(range(i+1, no), range(j+1, no)):
                count += 1

        for _ in combinations(range(no), 2):
            count += 1

        return [coeff] * count

    def mapping(self, profile, coeff):
        'Generate UCCGSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        value = iter(coeff)
        for i, j in combinations(range(no), 2):
            uccsd_fermion += singles(next(value), i, j)

        for i, j, k, l in combinations(range(no), 4):
            uccsd_fermion += doubles(next(value), i, j, k, l)
            uccsd_fermion += doubles(next(value), i+no, j+no, k+no, l+no)

        for i, j in product(range(no), repeat = 2):
            for k, l in product(range(i+1, no), range(j+1, no)):
                uccsd_fermion += doubles(next(value), i, j+no, k, l+no)

        for i, j in combinations(range(no), 2):
            uccsd_fermion += doubles(next(value), i, j+no, i+no, j)

        return JordanWignerMapper(uccsd_fermion)

class kUpCCGSD(UCCSD):
    'k-paired Unitary Coupled Cluster Generalized Singles and Doubles (kUpCCGSD) ansatz'
    def __init__(self, k=1):
        self.k = k

    def generate_coeff(self, profile, coeff=0.0):
        'Generate UCCGSD coefficients'
        no = profile.num_orb
        count = 0

        for _ in range(self.k):
            for _ in combinations(range(no), 2):
                count += 3

            for _ in combinations(range(no), 4):
                count += 2

            for i, j in product(range(no), repeat = 2):
                for _ in product(range(i+1, no), range(j+1, no)):
                    count += 1

        return [coeff] * count

    def mapping(self, profile, coeff):
        'Generate UCCGSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        value = iter(coeff)
        for _ in range(self.k):
            for i, j in combinations(range(no), 2):
                uccsd_fermion += singles(next(value), i, j)

            for i, j, k, l in combinations(range(no), 4):
                uccsd_fermion += doubles(next(value), i, j, k, l)
                uccsd_fermion += doubles(next(value), i+no, j+no, k+no, l+no)

            for i, j in product(range(no), repeat = 2):
                for k, l in product(range(i+1, no), range(j+1, no)):
                    uccsd_fermion += doubles(next(value), i, j+no, k, l+no)

            for i, j in combinations(range(no), 2):
                uccsd_fermion += doubles(next(value), i, i+no, j, j+no)

        return JordanWignerMapper(uccsd_fermion)
