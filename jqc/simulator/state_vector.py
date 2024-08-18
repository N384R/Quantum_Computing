from multiprocessing import Pool
import time
import numpy as np
from qiskit.quantum_info import Statevector, partial_trace
from qiskit import QuantumCircuit
from jqc.vqe.profile import Profile
from jqc.mapper.pauli_string import PauliString
from jqc.measure.angular_momentum import s_plus, s_minus, s_z

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

    def single_measure(self, args: tuple[np.ndarray, PauliString, complex]):
        'Measure the expectation value of a Pauli string'
        statevector, p_string, values = args
        if p_string.count_iden > 2:
            probability = self.get_rdm_trace(statevector, p_string)
        else:
            probability = statevector.conj().T @ p_string.matrix @ statevector
        expectation = float(probability.real) * values.real
        return expectation

    def get_overlap(self, state1, state2) -> float:
        'Get the square of the overlap between two states.'
        statevector1 = self.get_statevector(state1.circuit)
        statevector2 = self.get_statevector(state2.circuit)
        overlap_sq = abs(np.dot(statevector1.conj().T, statevector2))**2
        return overlap_sq

    @staticmethod
    def get_statevector(qc) -> np.ndarray:
        'Get the state vector of a quantum circuit.'
        statevector = Statevector(qc).data.reshape(-1, 1)
        real_part = np.where(abs(statevector.real) < 1e-15, 0, statevector)
        imag_part = np.where(abs(statevector.imag) < 1e-15, 0, statevector)
        return real_part + 1j * imag_part

    @staticmethod
    def get_rdm_trace(statevector, p_string):
        'Get the reduced density matrix of a quantum circuit.'
        reduce_idx = [idx for idx, pauli in p_string.items() if pauli.symbol == 'I']
        left_pauli = [pauli for pauli in p_string.values()   if pauli.symbol != 'I']
        reduced_density_matrix = partial_trace(statevector, reduce_idx).data
        reduced_operator = PauliString(left_pauli).matrix
        return np.trace(reduced_density_matrix @ reduced_operator)

    def measure_spin(self, profile: Profile) -> float:
        'Measure the spin of a quantum circuit.'
        s_x_and_s_y = s_plus(profile) * s_minus(profile) + s_minus(profile) * s_plus(profile)
        s_x_and_s_y_val = self.measure(profile.circuit, s_x_and_s_y, parallel=False)
        s_z_val = self.measure(profile.circuit, s_z(profile) * s_z(profile), parallel=False)
        return 0.5 * s_x_and_s_y_val + s_z_val
