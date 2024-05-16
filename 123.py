from itertools import permutations
import numpy as np
from pyscf import gto, scf, ao2mo

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.735', basis = 'sto-3g')
rhf = scf.RHF(mol)
EN = rhf.kernel()
HCORE = rhf.get_hcore()
D = rhf.make_rdm1()
S = rhf.get_ovlp()
C = rhf.mo_coeff
NUM = HCORE.shape[0]
HCORE = rhf.get_hcore()
TWO_ELEC = mol.intor('int2e')

HCORE_MO = C.T @ HCORE @ C
TWO_ELEC_MO = ao2mo.kernel(mol, C, TWO_ELEC, compact=False)
TWO_ELEC_MO = TWO_ELEC_MO.reshape((NUM, NUM, NUM, NUM)).transpose(0, 2, 1, 3)

FOCK_MO = np.zeros_like(HCORE)

for i in range(NUM):
    FOCK_MO[i, i] = HCORE_MO[i, i]

for i in range(NUM):
    for b in range(int(NUM/2)):
        print(f'sigma{i+1} =  h{i+1}{i+1} + 2J{i+1}{b+1} - K{i+1}{b+1}')
        FOCK_MO[i, i] += 2*TWO_ELEC_MO[i, b, i, b] - TWO_ELEC_MO[i, b, b, i]

print(FOCK_MO)
pyscf_mo_energy = rhf.mo_energy
calc_mo_energy = np.linalg.eigvalsh(FOCK_MO)

print(pyscf_mo_energy)
print(calc_mo_energy)

print(TWO_ELEC_MO)
