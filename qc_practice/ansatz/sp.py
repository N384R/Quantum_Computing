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
            for _ in range(no*2):
                count += 2

        for _ in range(no*2):
                count += 2

        return [coeff] * count

    @staticmethod
    def Agate(qc: QuantumCircuit, coeff, i):
        qc.cx(i+1, i)
        qc.ry(coeff, i)
        qc.cx(i, i+1)
        qc.ry(-coeff, i)
        qc.cx(i+1, i)

    def ansatz(self, qc, profile, coeff):
        'Generate SP ansatz'
        no = profile.num_orb
        value = iter(coeff)
        for _ in range(self.depth):
            for i in range(no):
                self.Agate(qc, next(value), i*2)
            
            for i in range(no):
                self.Agate(qc, next(value), i*2+1)


