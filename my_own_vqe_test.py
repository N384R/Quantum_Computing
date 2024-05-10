#%%
from qiskit import QuantumCircuit

ELECTRONS = 2
qc = QuantumCircuit(2*ELECTRONS)
for qubit in range(ELECTRONS):
    qc.x(qubit)

qc.draw('mpl')