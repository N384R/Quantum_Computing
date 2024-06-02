#%%
from pyscf import gto, scf

mol = gto.M(atom = 'Li 0 0 0; H 0 0 1.596', basis = 'sto-3g')
rhf = scf.RHF(mol)
e = rhf.kernel()
elec_energy = rhf.energy_elec()
print(elec_energy)

#%%
from pyscf import gto, scf

mol_o2 = gto.M(atom='O 0 0 0; O 0 0 1.2', spin=2) # (n+2 alpha, n beta) electrons
uhf_o2 = scf.UHF(mol_o2)
uhf_o2.kernel()
rohf_o2 = scf.ROHF(mol_o2)
rohf_o2.kernel()

# %%
from pyscf import gto, scf

mol = gto.M(
    atom = '''8  0  0.     0
              1  0  -0.757 0.587
              1  0  0.757  0.587''',
    basis = 'ccpvdz',
)

mf = scf.RHF(mol)
mf.conv_tol = 1e-1
mf.kernel()

#%%
from pyscf import gto, dft

mol = gto.M(atom = 'O 0 0 0; H 0 1 0; H 0 0 1', basis = 'ccpvdz')
rks_h2o = dft.RKS(mol)
rks_h2o.xc = 'b3lyp'
rks_h2o.kernel()

#%%
from pyscf import gto, scf, tdscf

mol_h2o = gto.M(atom = 'O 0 0 0; H 0 1 0; H 0 0 1', basis = 'ccpvdz')
rhf_h2o = scf.RHF(mol_h2o)
rhf_h2o.kernel()
tdhf_h2o = tdscf.TDHF(rhf_h2o)
tdhf_h2o.nstates = 6
tdhf_h2o.kernel()
weights, nto = tdhf_h2o.get_nto(state=2)

#%%
from pyscf import gto, scf, tdscf

mol_h2o = gto.M(atom = 'O 0 0 0; H 0 1 0; H 0 0 1', basis = 'ccpvdz')
rks_h2o = scf.RKS(mol_h2o)
rks_h2o.kernel()
tddft_h2o = tdscf.TDA(rks_h2o)
tddft_h2o.nstates = 6
tddft_h2o.kernel()
weights, nto = tddft_h2o.get_nto(state=2)

#%%

from pyscf import gto, scf
from pyscf.geomopt.geometric_solver import optimize

mol_h2o = gto.M(atom = 'O 0 0 0; H 0 1 0; H 0 0 1', basis = 'ccpvdz')
rhf_h2o = scf.RHF(mol_h2o)
opt_h2o = optimize(rhf_h2o)

# %%

from pyscf import gto, scf

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.75', basis = 'sto-3g')
mf = scf.dhf.UHF(mol)
mf.kernel()


# %%

from pyscf import gto, scf

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.75', basis = 'sto-3g')
mf = scf.UHF(mol)
mf.kernel()

#%%
from pyscf import gto, scf
from pyscf import cc

mol = gto.M(atom = 'H 0 0 0; H 0 0 2.5', basis = 'sto-3g')

mf = scf.UHF(mol)
mf.kernel()
ccsd = cc.CCSD(mf)
e = ccsd.kernel()
total_e = mf.e_tot + e[0]

nroots = 2

eomsf_singlet = ccsd.eomsf_ccsd(nroots=nroots)
print(eomsf_singlet[0])

print(f'Root 0: {total_e}')
for i in range(nroots):
    print(f'Root {i+1}: {total_e + eomsf_singlet[0][i]}')

# %%
from itertools import permutations
import numpy as np
from pyscf import gto, scf, ao2mo

mol = gto.M(atom = 'Li 0 0 0; H 0 0 0.735', basis = 'sto-3g')
rhf = scf.RHF(mol)
EN = rhf.kernel()
HCORE = rhf.get_hcore()
D = rhf.make_rdm1()
S = rhf.get_ovlp()
C = rhf.mo_coeff
electrons = mol.nelectron
NUM = HCORE.shape[0]
HCORE = rhf.get_hcore()
TWO_ELEC = mol.intor('int2e')

HCORE_MO = C.T @ HCORE @ C
TWO_ELEC_MO = ao2mo.kernel(mol, C, TWO_ELEC, compact=False)
TWO_ELEC_MO = TWO_ELEC_MO.reshape((NUM, NUM, NUM, NUM))

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

occ = rhf.mo_occ
print(occ)

#%%
import json
import numpy as np
from pyscf import gto, scf
from pyscf import cc

with open('data/H2_exact_PES.json', 'r', encoding='utf-8') as f:
    data = json.load(f)

for key, value in data.items():
    mol = gto.M(atom = f'H 0 0 0; H 0 0 {key}', basis = 'sto-3g', spin=2)
    mf = scf.UHF(mol)
    mf.kernel()
    ccsd = cc.CCSD(mf)
    e = ccsd.kernel()
    total_e = mf.e_tot + e[0]
    # total_e = mf.e_tot
    data[key].insert(1, total_e)
    print(f'\nBond Length: {key}')
    print('Energy:')
    print(data[key])

with open('data/H2_exact_PES.json', 'w', encoding='utf-8') as f:
    json.dump(data, f, indent=4)
