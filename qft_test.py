#%%
import numpy as np
from numpy import pi
from qiskit import QuantumCircuit
from qiskit.circuit.library import QFT
from qiskit.visualization import plot_histogram, plot_bloch_multivector
from qiskit_ibm_provider import IBMProvider

provider = IBMProvider()
simulator = provider.get_backend('ibmq_qasm_simulator')

num = 5

qc = QuantumCircuit(num+1, num)
qc.x(num)
for qubit in range(num):
    qc.h(qubit)

angle = 2*pi/3
for counting_qubit in range(num):
    for i in range(2**counting_qubit):
        qc.cp(angle, counting_qubit, num)

def qft_dagger(qc, n):
    for qubit in range(n//2):
        qc.swap(qubit, n-qubit-1)
    for j in range(n):
        for m in range(j):
            qc.cp(-pi/float(2**(j-m)), m, j)
        qc.h(j)

qc.barrier()
qft_dagger(qc, num)
qc.barrier()

for i in range(num):
    qc.measure(i, i)

qc.draw('mpl')


shots = 2048
job = simulator.run(qc, shots=shots)
result = job.result()
counts = result.get_counts(qc)

plot_histogram(counts)


# %%
