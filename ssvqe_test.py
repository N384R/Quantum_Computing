#%%
from pyscf import gto
from qc_practice import SSVQE

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.75', basis = 'sto-3g')
ssvqe = SSVQE(mol)
ssvqe.run()
