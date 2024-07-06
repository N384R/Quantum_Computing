from pyscf import gto
from qc_practice import VQE, SSVQE, VQD
from qc_practice.ansatz import UCCSD, kUpCCGSD
from qc_practice.simulator import StateVector

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.7', basis = '4-31g')
ssvqe = VQE(mol)
ssvqe.ansatz = UCCSD()
ssvqe.simulator = StateVector()
ssvqe.parallel = True

ssvqe.run()
