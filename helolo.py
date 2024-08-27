from functools import reduce
import numpy as np
from scipy.sparse import kron, csr_matrix
from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector, partial_trace
from jqc.mapper.pauli import Pauli

def matrix(pauli):
    return reduce(kron, [p.matrix for p in pauli])

I = Pauli('I')
X = Pauli('X')
Y = Pauli('Y')
Z = Pauli('Z')

p_string = (Y, Z, Z, Z, Z, Z, Z, Z, Y, I, I, X, Z, Z, Z, Z, X, X, X, I, I, I)
qc = QuantumCircuit(len(p_string))
qc.x(0)
qc.x(1)
qc.x(len(p_string)//2)
qc.x(len(p_string)//2 + 1)
state_vector = Statevector.from_instruction(qc)
reduce_idx = [idx for idx, p in enumerate(p_string) if p.symbol == 'I']
left_pauli = [p for p in p_string if p.symbol != 'I']

rdm = csr_matrix(partial_trace(state_vector, reduce_idx).data)
op = matrix(left_pauli)

mat = rdm @ op
trace = mat.trace()

print(f'size: {mat.toarray().nbytes / (1024**2)}MB')
