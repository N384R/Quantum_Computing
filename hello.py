#%%
import time
from pyscf import gto
from qc_practice import VQE
from qc_practice.ansatz import UCCSD, HEA

start = time.time()
mol = gto.M(atom = 'H 0 0 0; H 0 0 2.30', basis = 'sto-3g')
vqe = VQE(mol)
vqe.ansatz = UCCSD()
vqe.run()

print(f'VQE calculation took {time.time()-start:.2f} seconds.')

#%%

a = [1, 2, 3, 4]
b = sum(a)
print(b)