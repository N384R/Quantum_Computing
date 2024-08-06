#%%
import time
from pyscf import gto
from jqc import VQE, VQD
from jqc.ansatz import UCCSD, kUpCCGSD, UCCGSD

mol = gto.M(atom='Li 0 0 0; H 0 0 1.60', basis='sto-3g')
vqe = VQE(mol, UCCSD())
coeff = vqe.ansatz.generate_coeff(vqe.profile, coeff=1)
start = time.time()
vqe._batch(coeff)
print(f' Final batch time: {time.time()-start}')
# vqe.profile.circuit.draw()
# vqe.run()

#%%

from pyscf import gto
from jqc import VQE
from jqc.ansatz import kUpCCGSD, UCCGSD

mol = gto.M(atom='H 0 0 0; H 0 0 0.70', basis='4-31g')
vqe = VQE(mol, UCCGSD())
vqe.parallel = True

vqe.run()