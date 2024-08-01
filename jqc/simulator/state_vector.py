from multiprocessing import Pool
import numpy as np
from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit

class StateVector:
    'Class for running State Vector simulator.'

    def measure(self, qc: QuantumCircuit, operator, parallel: bool) -> float:
        'Measure the expectation value of a Hamiltonian'
        statevector = self.get_statevector(qc)
        tasks = [(statevector, p_string, values)
                 for p_string, values in operator.items()]
        if parallel:
            with Pool(4) as pool:
                energy = sum(pool.map(self.single_measure, tasks))
        else:
            energy = sum(self.single_measure(task) for task in tasks)
        return energy

    def single_measure(self, args: tuple[np.ndarray, dict, complex]):
        'Measure the expectation value of a Pauli string'
        statevector, p_string, values = args
        if all(p.symbol == 'I' for p in p_string.values()):
            return values.real
        probability = self.run_simulator(statevector, p_string)
        expectation = float(probability.real) * values.real
        return expectation

    def run_simulator(self, statevector, p_string) -> float:
        'Run the State Vector simulator.'
        probability = np.dot(statevector.conj().T, np.dot(p_string.matrix, statevector))
        return float(probability.real)

    def get_overlap(self, state1, state2) -> float:
        'Get the square of the overlap between two states.'
        statevector1 = self.get_statevector(state1.circuit)
        statevector2 = self.get_statevector(state2.circuit)
        overlap_sq = abs(np.dot(statevector1.conj().T, statevector2))**2
        return overlap_sq

    @staticmethod
    def get_statevector(qc) -> np.ndarray:
        'Get the state vector of a quantum circuit.'
        statevector: np.ndarray = np.asarray(Statevector(qc)).reshape(-1, 1)
        real_part = np.where(np.abs(statevector.real) < 1e-15, 0, statevector.real)
        imag_part = np.where(np.abs(statevector.imag) < 1e-15, 0, statevector.imag)
        statevector = real_part + 1j * imag_part
        return statevector
