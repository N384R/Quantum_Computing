from pyscf import gto
from qc_practice import VQD
from qc_practice.ansatz import UCCSD, UpCCGSD
from qc_practice.simulator import QASM

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.7', basis = 'sto-3g')
vqd = VQD(mol)
vqd.ansatz = UCCSD()
vqd.simulator = QASM()
vqd.parallel = True
vqd.run()