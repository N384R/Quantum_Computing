from qiskit_aer import AerProvider
from qiskit import QuantumCircuit
from qc_practice.ansatz import Ansatz

class QASM:
    'Class for running Quantum Assembly (QASM) simulator.'

    def __init__(self, shots = 10000):
        self.shots = shots
        self.backend = AerProvider().get_backend('qasm_simulator')

    def run_simulator(self, qc, p_string) -> float:
        'Run the QASM simulator.'
        for idx, p in p_string.items():
            if p.symbol == 'X':
                qc.h(idx)
            elif p.symbol == 'Y':
                qc.sdg(idx)
                qc.h(idx)
        for idx, p in p_string.items():
            if p.symbol == 'I':
                continue
            qc.measure(idx, idx)
        result = self.backend.run(qc, shots=self.shots).result().get_counts()
        counts = 0
        for key, value in result.items():  # type: ignore
            counts += (-1)**sum(int(k) for k in key) * value
        probability: float = counts / self.shots
        return probability

    def swap_test(self, state1, state2) -> float:
        'Swap test for measuring the overlap between two states.'
        no = state1.num_orb
        qc = circuit_swap_test(state1, state2)
        qc.measure(0, 0)
        result = self.backend.run(qc, shots=self.shots).result().get_counts()
        overlap_sq = abs(result.get('0'*(4*no + 1)) / self.shots * 2 - 1)
        return overlap_sq

def circuit_swap_test(state1, state2):
    'Swap test circuit for measuring the overlap between two states.'
    no = state1.num_orb
    qc: QuantumCircuit = QuantumCircuit(1, 1)
    for i, state in enumerate([state1, state2]):
        qc2 = state.circuit
        qc: QuantumCircuit = qc2.tensor(qc) # type: ignore
    qc.h(0)
    for i in range(1, no*2+1):
        qc.cswap(0, i, i+no*2)
    qc.h(0)
    return qc
