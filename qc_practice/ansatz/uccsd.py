'''
Unitary Coupled Cluster Single and Double (UCCSD) ansatz group

The class has the following methods:

generate_coeff: Generates UCCSD coefficients.
ansatz: Generates UCCSD ansatz circuit.
'''
from itertools import product
import numpy as np
import scipy.optimize as opt
from qc_practice.mapper.jordan_wigner import JordanWignerMapper

def sign_p(val):
    'Sign function'
    return '+' if val > 0 else '-'

def sign_m(val):
    'Sign function'
    return '-' if val > 0 else '+'

def boundary(coeff):
    'Boundary condition'
    return [(-np.pi, np.pi)] * len(coeff)

class UCCSD:
    'Unitary Coupled Cluster Single and Double (UCCSD) ansatz'

    @staticmethod
    def call_optimizer(func, coeff, method) -> opt.OptimizeResult:
        'Optimize the coefficients'
        return opt.minimize(func, coeff, method=method, bounds=boundary(coeff))

    def generate_coeff(self, profile, coeff=1e-3):
        'Generate UCCSD coefficients'
        no = profile.num_orb
        ne = profile.num_elec
        count = 0

        for _ in product(range(ne//2), range(ne//2, no)):
            count += 2
        for _ in product(range(ne//2), range(ne//2),
                        range(ne//2, no), range(ne//2, no)):
            count += 1
        return [coeff] * count

    def mapping(self, profile, coeff):
        'Generate UCCSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        ne = profile.num_elec
        value = iter(coeff)
        for i, j in product(range(ne//2), range(ne//2, no)):
            val = next(value)
            uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j}^ {i} ' + '\n'
            uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j} ' + '\n'

            val = next(value)
            uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j+no}^ {i+no} ' + '\n'
            uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i+no}^ {j+no} ' + '\n'

        for i, j, k, l in product(range(ne//2), range(ne//2),
                                  range(ne//2, no), range(ne//2, no)):
            val = next(value)
            uccsd_fermion += f'{sign_p(val)} {abs(val):f} {l}^ {k+no}^ {j+no} {i} ' + '\n'
            uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j+no}^ {k+no} {l} ' + '\n'
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

class UpCCSD(UCCSD):
    'paired Unitary Coupled Cluster Single and Double (UCCSD) ansatz'

    def generate_coeff(self, profile, coeff=1e-3):
        'Generate Spin Flip UCCSD coefficients'
        no = profile.num_orb
        ne = profile.num_elec
        count = 0
        for _ in product(range(ne//2), range(ne//2, no)):
            count += 2
        for _ in product(range(ne//2), range(ne//2),
                                  range(ne//2, no), range(ne//2, no)):
            count += 1
        for i in range(no-1):
            for _ in range(i+1, no):
                count += 1
        return [coeff] * count

    def mapping(self, profile, coeff):
        'Generate Spin Flip UCCSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        ne = profile.num_elec
        value = iter(coeff)
        for i in range(ne//2):
            for j in range(ne//2, no):
                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j}^ {i} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j} ' + '\n'

                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j+no}^ {i+no} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i+no}^ {j+no} ' + '\n'

        for i, j, k, l in product(range(ne//2), range(ne//2),
                                  range(ne//2, no), range(ne//2, no)):
            val = next(value)
            uccsd_fermion += f'{sign_p(val)} {abs(val):f} {l}^ {k+no}^ {j+no} {i} ' + '\n'
            uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j+no}^ {k+no} {l} ' + '\n'

        for i in range(no-1):
            for j in range(i+1, no):
                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j}^ {i+no}^ {j+no} {i} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j+no}^ {i+no} {j} ' + '\n'

        return JordanWignerMapper(uccsd_fermion)

class UCCGSD(UCCSD):
    'Unitary Coupled Cluster Generalized Single and Double (UCCGSD) ansatz'

    def generate_coeff(self, profile, coeff=1e-3):
        'Generate UCCGSD coefficients'
        no = profile.num_orb
        count = 0
        for i in range(no):
            for j in range(i+1, no):
                count += 2
        for i, j in product(range(no), range(no)):
            for _ in product(range(i+1, no), range(j+1, no)):
                count += 1
        return [coeff] * count

    def mapping(self, profile, coeff):
        'Generate UCCGSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        value = iter(coeff)
        for i in range(no):
            for j in range(i+1, no):
                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j}^ {i} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j} ' + '\n'

                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j+no}^ {i+no} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i+no}^ {j+no} ' + '\n'

        for i, j in product(range(no), range(no)):
            for k, l in product(range(i+1, no), range(j+1, no)):
                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {l}^ {k+no}^ {j+no} {i} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j+no}^ {k+no} {l} ' + '\n'

        return JordanWignerMapper(uccsd_fermion)

class UpCCGSD(UCCSD):
    'paired Unitary Coupled Cluster Generalized Single and Double (UpCCGSD) ansatz'

    def generate_coeff(self, profile, coeff=1e-3):
        'Generate UCCGSD coefficients'
        no = profile.num_orb
        count = 0

        for i in range(no):
            for j in range(i+1, no):
                count += 2

        for i, j in product(range(no), range(no)):
            for _ in product(range(i+1, no), range(j+1, no)):
                count += 1

        for i in range(no-1):
            for _ in range(i+1, no):
                count += 1

        return [coeff] * count

    def mapping(self, profile, coeff):
        'Generate UCCGSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        value = iter(coeff)

        for i in range(no):
            for j in range(i+1, no):
                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j}^ {i} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j} ' + '\n'

                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j+no}^ {i+no} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i+no}^ {j+no} ' + '\n'

        for i, j in product(range(no), range(no)):
            for k, l in product(range(i+1, no), range(j+1, no)):
                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {l}^ {k+no}^ {j+no} {i} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j+no}^ {k+no} {l} ' + '\n'

        for i in range(no):
            for j in range(i+1, no):
                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j}^ {i+no}^ {j+no} {i} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j+no}^ {i+no} {j} ' + '\n'

        return JordanWignerMapper(uccsd_fermion)

class kUpCCGSD(UCCSD):
    'k-paired Unitary Coupled Cluster Generalized Single and Double (kUpCCGSD) ansatz'
    def __init__(self, k=1):
        self.k = k

    def generate_coeff(self, profile, coeff=1e-3):
        'Generate UCCGSD coefficients'
        no = profile.num_orb
        ne = profile.num_elec
        count = 0

        for i in range(no):
            for j in range(i+1, no):
                count += 2

        for i, j in product(range(ne//2 + self.k), range(ne//2 + self.k)):
            for _ in product(range(i+1, no), range(j+1, no)):
                count += 1

        for i in range(ne//2 + self.k):
            for _ in range(i+1, no):
                count += 1

        return [coeff] * count

    def mapping(self, profile, coeff):
        'Generate UCCGSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        ne = profile.num_elec
        value = iter(coeff)
        for i in range(no):
            for j in range(i+1, no):
                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j}^ {i} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j} ' + '\n'

                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j+no}^ {i+no} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i+no}^ {j+no} ' + '\n'

        for i, j in product(range(ne//2 + self.k), range(ne//2 + self.k)):
            for k, l in product(range(i+1, no), range(j+1, no)):
                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {l}^ {k+no}^ {j+no} {i} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j+no}^ {k+no} {l} ' + '\n'

        for i in range(ne//2 + self.k):
            for j in range(i+1, no):
                val = next(value)
                uccsd_fermion += f'{sign_p(val)} {abs(val):f} {j}^ {i+no}^ {j+no} {i} ' + '\n'
                uccsd_fermion += f'{sign_m(val)} {abs(val):f} {i}^ {j+no}^ {i+no} {j} ' + '\n'

        return JordanWignerMapper(uccsd_fermion)
