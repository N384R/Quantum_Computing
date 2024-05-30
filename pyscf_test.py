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