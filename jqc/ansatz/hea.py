'''
Hardware Efficient Ansatz (HEA) for VQE

The class has the following methods:

generate_coeff: Generates HEA coefficients.
ansatz: Generates HEA ansatz circuit.
'''
import numpy as np
import scipy.optimize as opt


def boundary(coeff):
    'Boundary condition'
    return [(-2*np.pi, 2*np.pi)] * len(coeff)


class HEA:
    'Hardware Efficient Ansatz (HEA)'

    def __init__(self, depth=1):
        self.depth = depth

    @staticmethod
    def call_optimizer(func, coeff, method):
        'Optimize the coefficients'
        return opt.minimize(func, coeff, method=method, bounds=boundary(coeff))

    def generate_coeff(self, profile, coeff=0.0):
        'Generate HEA coefficients'
        no = profile.num_orb
        count = 0
        for _ in range(self.depth):
            for _ in range(no*2):
                count += 2

        for _ in range(no*2):
            count += 2

        return [coeff] * count

    def ansatz(self, qc, profile, coeff):
        'Generate HEA ansatz'
        no = profile.num_orb
        value = iter(coeff)
        for _ in range(self.depth):
            for i in range(no*2):
                qc.ry(next(value), i)
                qc.rz(next(value), i)

            for i in range(no*2-1):
                qc.cx(i, i+1)

        for i in range(no*2):
            qc.ry(next(value), i)
            qc.rz(next(value), i)
