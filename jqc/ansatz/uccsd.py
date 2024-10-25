'''
Unitary Coupled Cluster Singles and Doubles (UCCSD) ansatz group

The class has the following methods:

generate_coeff: Generates UCCSD coefficients.
ansatz: Generates UCCSD ansatz circuit.
'''
from itertools import product, combinations
import numpy as np
import scipy.optimize as opt
from jqc.mapper.fermion import FermionicOp


def singles(val, i, j):
    'Singles operator'
    return (FermionicOp(val, f'{j}^ {i}') -
            FermionicOp(val, f'{i}^ {j}'))


def doubles(val, i, j, k, l):
    'Doubles operator'
    return (FermionicOp(val, f'{k}^ {l}^ {j} {i}') -
            FermionicOp(val, f'{i}^ {j}^ {l} {k}'))


def boundary(coeff):
    'Boundary condition'
    return [(-2*np.pi, 2*np.pi)] * len(coeff)


class UCCSD:
    'Unitary Coupled Cluster Singles and Doubles (UCCSD) ansatz'

    @staticmethod
    def call_optimizer(func, coeff, method) -> opt.OptimizeResult:
        'Optimize the coefficients'
        if method in ('BFGS', 'CG'):
            return opt.minimize(func, coeff, method=method)
        return opt.minimize(func, coeff, method=method, bounds=boundary(coeff))

    def generate_coeff(self, profile, coeff=0.0):
        'Generate coefficients'
        count = self._term_count(profile)
        return [coeff] * count

    def _term_count(self, profile):
        'Count the number of terms in the UCCSD ansatz'
        no = profile.num_orb
        ne = profile.num_elec // 2
        nsingles = len(list(self._singles_idx(ne, no)))
        ndoubles = len(list(self._doubles_idx(ne, no)))
        return nsingles + ndoubles

    def _singles_idx(self, ne, no):
        'Generate singles indices'
        for i in range(ne):
            for j in range(ne, no):
                yield i, j
                yield i+no, j+no

    def _doubles_idx(self, ne, no):
        'Generate doubles indices'
        for i, j in product(range(ne), repeat=2):
            for k, l in product(range(ne, no), repeat=2):
                yield i, j+no, k, l+no

    def mapping(self, profile, coeff):
        'Generate ansatz Pauli strings'
        fermionic_op = FermionicOp()
        no = profile.num_orb
        ne = profile.num_elec // 2
        value = iter(coeff)
        for p, q in self._singles_idx(ne, no):
            fermionic_op += singles(next(value), p, q)

        for p, q, r, s, in self._doubles_idx(ne, no):
            fermionic_op += doubles(next(value), p, q, r, s)
        return fermionic_op.jordan_wigner

    def ansatz(self, qc, profile, coeff):
        'Generate UCCSD ansatz circuit'
        for obj, val in self.mapping(profile, coeff).items():
            chk = []
            for idx, po in enumerate(obj):
                if po.symbol == 'X':
                    qc.h(idx)
                elif po.symbol == 'Y':
                    qc.sdg(idx)
                    qc.h(idx)
                if po.symbol != 'I':
                    chk.append(idx)

            for i, j in zip(chk, chk[1:]):
                qc.cx(i, j)

            qc.rz(val.imag, max(chk))

            for i, j in zip(chk[:-1][::-1], chk[1:][::-1]):
                qc.cx(i, j)

            for idx, po in enumerate(obj):
                if po.symbol == 'X':
                    qc.h(idx)
                elif po.symbol == 'Y':
                    qc.h(idx)
                    qc.s(idx)


class fUCCSD(UCCSD):
    ' spin-flip Unitary Coupled Cluster Singles and Doubles (fUCCSD) ansatz'

    def _singles_idx(self, ne, no):
        for i in range(ne):
            for j in range(ne, no):
                yield i, j
                yield i+no, j+no
                yield i, j+no
                yield i+no, j

    def _doubles_idx(self, ne, no):
        yield from super()._doubles_idx(ne, no)
        for i, j in combinations(range(no), 2):
            yield i, j+no, i+no, j


class UCCGSD(UCCSD):
    'Unitary Coupled Cluster Generalized Singles and Doubles (UCCGSD) ansatz'

    def _singles_idx(self, ne, no):
        'Generate singles indices'
        for i, j in combinations(range(no), 2):
            yield i, j
            yield i+no, j+no

    def _doubles_idx(self, ne, no):
        'Generate doubles indices'
        for i, j, k, l in combinations(range(no), 4):
            yield i, j, k, l
            yield i+no, j+no, k+no, l+no

        for i, j in product(range(no), repeat=2):
            for k, l in product(range(i+1, no), range(j+1, no)):
                yield i, j+no, k, l+no


class eUCCGSD(UCCGSD):
    'entangled Unitary Coupled Cluster Generalized Singles and Doubles (eUCCGSD) ansatz'

    def _doubles_idx(self, ne, no):
        yield from super()._doubles_idx(ne, no)
        for i, j in combinations(range(no), 2):
            yield j, i, i, j


class kUpCCGSD(UCCGSD):
    'k-paired Unitary Coupled Cluster Generalized Singles and Doubles (kUpCCGSD) ansatz'

    def __init__(self, k=1):
        self.k = k

    def _singles_idx(self, ne, no):
        'Generate singles indices'
        for _ in range(self.k):
            yield from super()._singles_idx(ne, no)

    def _doubles_idx(self, ne, no):
        'Generate doubles indices'
        for _ in range(self.k):
            yield from super()._doubles_idx(ne, no)
