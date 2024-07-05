from pyscf import gto
from qc_practice.ssvqe import SSVQE
from qc_practice.ansatz import UCCSD, UpCCGSD

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.7', basis = 'sto-3g')
ssvqe = SSVQE(mol)
ssvqe.ansatz = UCCSD()

ssvqe.run()