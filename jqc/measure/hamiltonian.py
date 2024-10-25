from itertools import product
import numpy as np
from pyscf import ao2mo, scf
from jqc.mapper.fermion import FermionicOp


def pyscf_luncher(mol):
    '''Run PySCF to get the Hamiltonian information.'''

    mf = scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()

    c = np.asarray(mf.mo_coeff)
    hcore_ao = mf.get_hcore()
    two_elec_ao = mol.intor('int2e')
    two_elec_mo = np.asarray(ao2mo.kernel(mol, c, two_elec_ao, compact=False)
                             ).reshape((mol.nao, mol.nao, mol.nao, mol.nao))
    dip_ao = mol.intor('int1e_r', comp=3)
    return {
        'num_orb': mol.nao,
        'num_elec': mol.nelectron,
        'hcore_mo': c.T @ hcore_ao @ c,
        'two_elec_mo': two_elec_mo,
        'dipole_mo': np.einsum('xuv,up,vq->xpq', dip_ao, c, c),
        'energy_nuc': mol.energy_nuc(),
        'energy_elec': mf.energy_elec()[0]
    }


def hamiltonian(profile):
    '''Generate the Hamiltonian in the second quantization form.'''
    second_q = FermionicOp()
    n = profile.num_orb
    for i, j in product(range(n), repeat=2):
        if abs(coeff := profile.qm['hcore'][i, j]) < 1e-10:
            continue
        second_q += FermionicOp(coeff, f'{j}^ {i}') + \
            FermionicOp(coeff, f'{j+n}^ {i+n}')
    for i, j, k, l in product(range(n), repeat=4):
        if abs(coeff := profile.qm['two_elec'][i, j, k, l]/2) < 1e-10:
            continue
        second_q += (FermionicOp(coeff, f'{k}^ {j}^ {i} {l}') +
                     FermionicOp(coeff, f'{k}^ {j+n}^ {i+n} {l}') +
                     FermionicOp(coeff, f'{k+n}^ {j}^ {i} {l+n}') +
                     FermionicOp(coeff, f'{k+n}^ {j+n}^ {i+n} {l+n}'))
    return second_q.jordan_wigner
