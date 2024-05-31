#%%
from pyscf import gto
from qc_practice import VQE

mol = gto.M(atom = 'Li 0 0 0; H 0 0 1.596', basis = 'sto-3g')
vqe = VQE(mol)
vqe.run()

#%%
from pyscf import gto
from qc_practice import VQE

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.75', basis = 'sto-3g')
vqe = VQE(mol)
vqe.run()

#%%

from pyscf import gto
from qc_practice import SSVQE

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.74', basis = 'sto-3g')
ssvqe = SSVQE(mol)
ssvqe.verbose = 0
ssvqe.run()
