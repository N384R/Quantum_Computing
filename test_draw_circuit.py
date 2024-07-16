#%%
from pyscf import gto
from qc_practice import VQE, VQD
from qc_practice.ansatz import UCCSD, kUpCCGSD, UCCGSD

mol = gto.M(atom='H 0 0 0; H 0 0 0.70', basis='4-31g')
vqe = VQE(mol, UCCGSD())
coeff = vqe.ansatz.generate_coeff(vqe.profile, coeff=1e-5)
vqe._batch(coeff)
vqe.profile.circuit.draw()
# vqe.run()