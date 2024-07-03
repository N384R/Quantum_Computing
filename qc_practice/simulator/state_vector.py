from qiskit_aer import AerProvider
from qiskit import QuantumCircuit
from qc_practice.ansatz import Ansatz
from qc_practice.simulator.qasm import circuit_swap_test

class StateVector:  ####need to build correctly
    'Class for running State Vector simulator.'

    def __init__(self):
        self.backend = AerProvider().get_backend('statevector_simulator')

    def run_simulator(self, qc, p_string) -> float:
        'Run the State Vector simulator.'
        result = self.backend.run(qc).result().get_statevector()
        return result

    def swap_test(self, state1, state2) -> float:
        'Swap test for measuring the overlap between two states.'
        qc = circuit_swap_test(state1, state2)
        result = self.backend.run(qc).result().get_statevector()
        overlap_sq = abs(result[0] * result[0].conj())
        return overlap_sq

