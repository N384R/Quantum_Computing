from pyscf import gto
from jqc import VQE, VQD, SSVQE
from jqc.ansatz import fUCCSD, SP

mol = gto.M(atom='H 0 0 0; H 0 0 0.70', basis='4-31g')
ssvqe = SSVQE(mol, SP(depth=7))
ssvqe.active_space = [1, 1]

ssvqe.run()
