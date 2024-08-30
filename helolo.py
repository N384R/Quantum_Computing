'병렬 있음'

from pyscf import gto
from jqc import VQE
from jqc.ansatz import UCCSD

mol = gto.M(atom='Li 0 0 0; H 0 0 1.60', basis='sto-3g')
vqe = VQE(mol, UCCSD())
vqe.parallel = True

vqe.run()