#%%
from pyscf import gto
from qc_practice import SSVQE
from qc_practice.ansatz import SP
from qiskit import QuantumCircuit
from qc_practice.profile import Profile


sp = SP()
profile = Profile()
profile.num_orb = 
qc = QuantumCircuit(profile.num_orb*2)
coeff = sp.generate_coeff(profile)
sp.ansatz(qc, profile, coeff)
qc.draw('mpl')