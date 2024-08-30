from joblib import Parallel, delayed
import numpy as np
from numba import jit
from qiskit.quantum_info import Statevector
from qiskit import QuantumCircuit
from jqc.mapper.pauli import Pauli

class StateVector:
    'Class for running State Vector simulator.'

    def measure(self, qc: QuantumCircuit, operator, parallel = False) -> float:
        'Measure the expectation value of an Operator'
        sv = get_statevector(qc)
        tasks = ((sv, p_string, values) for p_string, values in operator.items())
        if parallel:
            quantities = Parallel(n_jobs=-1)(delayed(self.single_measure)(task) for task in tasks)
            quantity = sum(quantities) # type: ignore
        else:
            quantity = sum(self.single_measure(task) for task in tasks)
        return quantity

    @staticmethod
    def single_measure(args: tuple) -> float:
        'Measure the expectation value of a Pauli string'
        statevector, p_string, values = args
        if np.all(p_string == Pauli.I):
            return values.real
        statevector2 = mult_operator(statevector, pauli_to_int(p_string))
        probability = np.dot(statevector.conj().T, statevector2)
        expectation = float(probability.real) * values.real
        return expectation

    def get_overlap(self, state1, state2) -> float:
        'Get the square of the overlap between two states.'
        statevector1 = get_statevector(state1.circuit)
        statevector2 = get_statevector(state2.circuit)
        overlap_sq = abs(np.dot(statevector1.conj().T, statevector2))**2
        return overlap_sq

def get_statevector(qc) -> np.ndarray:
    'Get the state vector of a quantum circuit.'
    statevector = Statevector(qc).data
    real_part = np.where(abs(statevector.real) < 1e-12, 0, statevector.real)
    imag_part = np.where(abs(statevector.imag) < 1e-12, 0, statevector.imag)
    return real_part + 1j * imag_part

def pauli_to_int(pauli_string):
    'Convert a Pauli string to an integer.'
    result = np.zeros(len(pauli_string), dtype=int)
    for i, p in enumerate(pauli_string):
        if p == Pauli.I:
            result[i] = 0
        elif p == Pauli.X:
            result[i] = 1
        elif p == Pauli.Y:
            result[i] = 2
        elif p == Pauli.Z:
            result[i] = 3
    return result

@jit(nopython=True)
def mult_operator(statevector, operator):
    'Multiply an operator to a statevector.'
    indices = np.arange(len(statevector))
    result = np.copy(statevector)
    for i, op in enumerate(operator):
        if op == 0:
            continue
        if op == 1:
            new_indices = indices ^ (1 << i)
            result[new_indices] = result[indices]
        elif op == 2:
            new_indices = indices ^ (1 << i)
            phase_factors = np.where((indices >> i) & 1, 1j, -1j)
            result[new_indices] = phase_factors * result[indices]
        elif op == 3:
            phase_factors = np.where((indices >> i) & 1, -1, 1)
            result *= phase_factors
    return result
