import numpy as np
from qiskit.quantum_info import Statevector

class StateVector:  ####need to build correctly
    'Class for running State Vector simulator.'
    def run_simulator(self, qc, p_string) -> float:
        'Run the State Vector simulator.'
        statevector = self.get_statevector(qc)
        probability = np.dot(statevector.conj().T, np.dot(p_string.matrix, statevector))
        return float(probability.real)

    def get_overlap(self, state1, state2) -> float:
        'Get the square of the overlap between two states.'
        statevector1 = self.get_statevector(state1.circuit)
        statevector2 = self.get_statevector(state2.circuit)
        overlap_sq = abs(np.dot(statevector1.conj().T, statevector2))**2
        return overlap_sq

    @staticmethod
    def get_statevector(qc):
        'Get the state vector of a quantum circuit.'
        statevector: np.ndarray = np.asarray(Statevector(qc)).reshape(-1, 1)
        real_part = np.where(np.abs(statevector.real) < 1e-15, 0, statevector.real)
        imag_part = np.where(np.abs(statevector.imag) < 1e-15, 0, statevector.imag)
        statevector = real_part + 1j * imag_part
        return statevector
