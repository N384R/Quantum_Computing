import numpy as np
from qiskit import QuantumCircuit

class SP:
    'Symmetry Preserving (SP) ansatz'
    def __init__(self, depth=1):
        self.depth = depth

    def generate_coeff(self, profile, coeff=1e-5):
        'Generate SP ansatz coefficients'
        no = profile.num_orb
        count = 0
        for _ in range(self.depth):
            for _ in range(0, no*2-1, 2):
                count += 2

            for _ in range(1, no*2-1, 2):
                count += 2

        for _ in range(0, no*2-1, 2):
            count += 2

        return [coeff] * count

    @staticmethod
    def Agate(qc: QuantumCircuit, val1, val2, i):
        'Generate Agate'
        qc.cx(i+1, i)
        qc.ry(-val2 - np.pi/2, i+1)
        qc.rz(-val1 - np.pi, i+1)
        qc.cx(i, i+1)
        qc.rz(val1 + np.pi, i+1)
        qc.ry(val2 + np.pi/2, i+1)
        qc.cx(i+1, i)

    def ansatz(self, qc, profile, coeff):
        'Generate SP ansatz'
        no = profile.num_orb
        value = iter(coeff)
        for _ in range(self.depth):
            for i in range(0, no*2-1, 2):
                self.Agate(qc, next(value), next(value), i)

            for i in range(1, no*2-1, 2):
                self.Agate(qc, next(value), next(value), i)

        for i in range(0, no*2-1, 2):
            self.Agate(qc, next(value), next(value), i)

class RSP:
    'Real-valued Symmetry Preserving (RSP) ansatz'
    def __init__(self, depth=1):
        self.depth = depth

    def generate_coeff(self, profile, coeff=1e-5):
        'Generate RSP ansatz coefficients'
        no = profile.num_orb
        count = 0
        for _ in range(self.depth):
            for _ in range(0, no*2-1, 2):
                count += 1

            for _ in range(1, no*2-1, 2):
                count += 1

        for _ in range(0, no*2-1, 2):
            count += 1

        return [coeff] * count

    @staticmethod
    def gate(qc: QuantumCircuit, val, i):
        'Generate gate'
        qc.cx(i+1, i)
        qc.ry(-val - np.pi/2, i+1)
        qc.cx(i, i+1)
        qc.ry(val + np.pi/2, i+1)
        qc.cx(i+1, i)

    def ansatz(self, qc, profile, coeff):
        'Generate RSP ansatz'
        no = profile.num_orb
        value = iter(coeff)
        for _ in range(self.depth):
            for i in range(0, no*2-1, 2):
                self.gate(qc, next(value), i)

            for i in range(1, no*2-1, 2):
                self.gate(qc, next(value), i)

        for i in range(0, no*2-1, 2):
            self.gate(qc, next(value), i)
