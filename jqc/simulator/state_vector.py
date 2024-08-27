from multiprocessing import Pool
import time
import numpy as np
from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit
from jqc.vqe.profile import Profile
from jqc.mapper.pauli import Pauli
from jqc.measure.angular_momentum import s_plus, s_minus, s_z


class StateVector:
    'Class for running State Vector simulator.'

    def measure(self, qc: QuantumCircuit, operator, parallel = False) -> float:
        'Measure the expectation value of an Operator'
        sv = get_statevector(qc)
        tasks = [(sv, p_string, values) for p_string, values in operator.items()]

        if parallel:
            with Pool(8) as pool:
                chunk_size = len(tasks) // pool._processes # type: ignore # pylint: disable=protected-access
                energies = pool.map(self.single_measure, tasks, chunksize=chunk_size)
            energy = sum(energies)
        else:
            energy = sum(self.single_measure(task) for task in tasks)
        return energy

    @staticmethod
    def single_measure(args: tuple[np.ndarray, tuple, complex]):
        'Measure the expectation value of a Pauli string'
        statevector, p_string, values = args
        if all(p == Pauli.I for p in p_string):
            return values.real
        statevector2 = mult_operator(statevector, p_string)
        probability = np.dot(statevector.conj().T, statevector2)
        expectation = float(probability.real) * values.real
        return expectation

    def get_overlap(self, state1, state2) -> float:
        'Get the square of the overlap between two states.'
        statevector1 = get_statevector(state1.circuit)
        statevector2 = get_statevector(state2.circuit)
        overlap_sq = abs(np.dot(statevector1.conj().T, statevector2))**2
        return overlap_sq

    def measure_spin(self, profile: Profile) -> float:
        'Measure the spin of a quantum circuit.'
        s_x_and_s_y = s_plus(profile) * s_minus(profile) + s_minus(profile) * s_plus(profile)
        s_x_and_s_y_val = self.measure(profile.circuit, s_x_and_s_y, parallel=False)
        s_z_val = self.measure(profile.circuit, s_z(profile) * s_z(profile), parallel=False)
        return 0.5 * s_x_and_s_y_val + s_z_val

def get_statevector(qc) -> np.ndarray:
    'Get the state vector of a quantum circuit.'
    statevector = Statevector(qc).data
    real_part = np.where(abs(statevector.real) < 1e-15, 0, statevector.real)
    imag_part = np.where(abs(statevector.imag) < 1e-15, 0, statevector.imag)
    # if np.all(imag_part == 0):
    #     return real_part
    return real_part + 1j * imag_part

def mult_operator(statevector, operator):
    'Multiply an operator to a statevector.'
    indices = np.arange(len(statevector))
    result = np.copy(statevector)
    for i, op in enumerate(operator):
        if op == Pauli.I:
            continue
        if op == Pauli.X:
            new_indices = indices ^ (1 << i)
            result[new_indices] = result[indices]
        elif op == Pauli.Y:
            new_indices = indices ^ (1 << i)
            result = result.astype(complex)
            phase_factors = np.where((indices >> i) & 1, 1j, -1j)
            result[new_indices] = phase_factors * result[indices]
        elif op == Pauli.Z:
            result = result.copy()
            phase_factors = np.where((indices >> i) & 1, -1, 1)
            result *= phase_factors
    return result
