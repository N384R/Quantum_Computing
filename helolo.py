from pyscf import gto
from jqc import VQE, VQD, SSVQE
from jqc.ansatz import eUCCSD

mol = gto.M(atom='H 0 0 0; H 0 0 0.70', basis='sto-3g')
ssvqe = SSVQE(mol, eUCCSD())
ssvqe.active_space = [1, 1]

ssvqe.run()