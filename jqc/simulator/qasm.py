from joblib import Parallel, delayed
from qiskit_aer import AerProvider
from qiskit import QuantumCircuit
from jqc.mapper.pauli import Pauli

class QASM:
    'Class for running Quantum Assembly (QASM) simulator.'

    def __init__(self, shots = 10000):
        self.shots = shots
        self.parallel = False
        self.backend = AerProvider().get_backend('qasm_simulator')

    def measure(self, qc1, operator, qc2 = None) -> float:
        'Measure the expectation value of a Hamiltonian'
        tasks = [(qc1, qc2, p_string, values) for p_string, values in operator.items()]
        if self.parallel:
            quantities = Parallel(n_jobs=-1)(delayed(self.single_measure)(task) for task in tasks)
            quantity = sum(quantities) # type: ignore
        else:
            quantity = sum(self.single_measure(task) for task in tasks)
        return quantity

    def single_measure(self, args: tuple) -> float:
        'Measure the expectation value of a Pauli string'
        qc1, qc2, p_string, values = args
        if all(p == Pauli.I for p in p_string.values()):
            return values.real
        probability = self.run_simulator(qc1, p_string)
        expectation = float(probability.real) * values.real
        return expectation

    def run_simulator(self, qc, p_string) -> float:
        'Run the QASM simulator.'
        for idx, p in p_string.items():
            if p == Pauli.X:
                qc.h(idx)
            elif p == Pauli.Y:
                qc.sdg(idx)
                qc.h(idx)
        for idx, p in p_string.items():
            if p == Pauli.I:
                continue
            qc.measure(idx, idx)
        result = self.backend.run(qc, shots=self.shots).result().get_counts()
        counts = 0
        for key, value in result.items():  # type: ignore
            counts += (-1)**sum(int(k) for k in key) * value
        probability: float = counts / self.shots
        return probability

    def get_overlap(self, state1, state2) -> float:
        'Get the square of the overlap between two states.'
        return self.swap_test(state1, state2)

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
