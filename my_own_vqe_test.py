 #%%
from pyscf import gto
from qc_practice import VQE

mol = gto.M(atom = 'Li 0 0 0; H 0 0 1.596', basis = 'sto-3g')
vqe = VQE(mol)
vqe.run()

#%%
from pyscf import gto
from qc_practice import VQE
from qc_practice.ansatz import UCCSD

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.75', basis = 'sto-3g')
vqe = VQE(mol, UCCSD())
vqe.run()

#%%

from pyscf import gto
from qc_practice import SSVQE
from qc_practice.ansatz import UCCSD, SpinFlipUCCSD

mol = gto.M(atom = 'H 0 0 0; H 0 0 1.90', basis = 'sto3-g')
ssvqe = SSVQE(mol)
ssvqe.weights = [1, 0.2, 0.2, 0.2, 0.05, 0.01]
ssvqe.ansatz = SpinFlipUCCSD()
ssvqe.run()

#%%

from pyscf import gto
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import SSVQE
from qc_practice.ansatz import SpinFlipUCCSD


mol = gto.M(atom = 'H 0 0 0; H 0 0 0.74', basis = 'sto-3g')

ssvqe = SSVQE(mol)
ssvqe.ansatz = SpinFlipUCCSD()
ssvqe.weights = [1, 0.5, 0.5, 0.5, 0.1, 0.01]
ssvqe.run()

for profile in ssvqe.profile:
    qc = profile.circuit
    energy = profile.energy_total()

    print('Spin measurement')
    qc.barrier()
    for k in range(4):
        qc.measure(k, k)

    qc.draw('mpl')

    # shots = 100000
    # backend = AerProvider().get_backend('qasm_simulator')
    # result = backend.run(qc, shots=shots).result().get_counts()
    # print(result)

    # spin = 0
    # for key, value in result.items():
    #     for i, orb in enumerate(reversed([*key])):
    #         if i < 2:
    #             spin += value * int(orb)
    #         else:
    #             spin -= value * int(orb)

    # print(f'multiplicity: {abs(spin/shots/2)}')

#%%
import numpy as np
from pyscf import gto, scf
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import SSVQE
from qc_practice.ansatz import SpinFlipUCCSD
from qc_practice.profile import Profile

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.74', basis = 'sto-3g')
rhf = scf.RHF(mol)
profile = Profile()
profile.num_orb = rhf.get_hcore().shape[0]

ssvqe = SSVQE(mol)
ssvqe.ansatz = SpinFlipUCCSD()
coeff = [0.01201245,  0.01782515, -0.23119888,  1.79793851]
qc = QuantumCircuit(4, 4)
qc.x(0)
qc.x(3)

ssvqe._circuit(qc, SpinFlipUCCSD().ansatz(profile, coeff))
qc.barrier()
for k in range(4):
    qc.measure(k, k)

qc.draw('mpl')
# shots = 100000
# backend = AerProvider().get_backend('qasm_simulator')
# result = backend.run(qc, shots=shots).result().get_counts()
# print(result)

# spin = 0
# for key, value in result.items():
#     for i, orb in enumerate(reversed([*key])):
#         if i < 2:
#             spin += value * int(orb)
#         else:
#             spin -= value * int(orb)

# print(f'multiplicity: {abs(spin/shots/2)}')