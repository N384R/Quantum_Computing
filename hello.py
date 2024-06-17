#%%
from pyscf import gto
from qc_practice import SSVQE
from qc_practice.ansatz import SP, RSP
from qiskit import QuantumCircuit
from qc_practice.profile import Profile

sp = SP(depth=3)
profile = Profile()
profile.num_orb = 2
qc = QuantumCircuit(profile.num_orb*2)
coeff = sp.generate_coeff(profile)
sp.ansatz(qc, profile, coeff)
qc.draw('mpl')

#%%

from qiskit import QuantumCircuit
from qc_practice.profile import Profile

profile = Profile()
profile.num_orb = 2
qc = QuantumCircuit(profile.num_orb*2, profile.num_orb*2)
qc_ancila = QuantumCircuit(1, 1)
qc_ancila.h(0)
qc_2 = qc.tensor(qc)
qc_2.tensor(qc_ancila, inplace=True)

qc_2.draw('mpl')

#%%

from pyscf import gto
from qc_practice import VQD, SSVQE
from qc_practice.ansatz import SP

mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='sto-3g')
vqd = VQD(mol)
vqd.ansatz = SP(depth=1)
vqd.run()
