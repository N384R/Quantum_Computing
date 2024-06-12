from itertools import product
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
            for i in range(0, no*2-1, 2):
                count += 2

            for i in range(1, no*2-1, 2):
                count += 2

        return [coeff] * count

    @staticmethod
    def Agate(qc: QuantumCircuit, coeff1, coeff2, i):
        'Generate Agate'
        qc.cx(i+1, i)
        qc.ry(coeff2, i+1)
        qc.rz(-coeff1, i+1)
        qc.cx(i, i+1)
        qc.rz(coeff1, i+1)
        qc.ry(coeff2, i+1)
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

