#%%
from pyscf import gto
from qc_practice import VQE, VQD
from qc_practice.ansatz import UCCSD, UpCCSD

mol = gto.M(atom='H 0 0 0; H 0 0 0.70', basis='4-31g')
vqe = VQE(mol, UCCSD())
vqe.parallel = True
coeff = vqe.ansatz.generate_coeff(vqe.profile)
vqe._batch(coeff)
vqe.profile.circuit.draw('mpl')
# vqe.run()
