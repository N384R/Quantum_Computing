#%%
from pyscf import gto
from qc_practice import VQE

mol = gto.M(atom = 'Li 0 0 0; H 0 0 1.596', basis = 'sto-3g')
vqe = VQE(mol)
vqe.run()
