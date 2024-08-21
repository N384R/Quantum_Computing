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
    hcore = mf.get_hcore()
    two_elec = mol.intor('int2e')
    result = {
        'hcore_mo': c.T @ hcore @ c,
        'num_elec': mol.nelectron,
        'num_orb' : hcore.shape[0],
        'two_elec_mo': np.asarray(ao2mo.kernel(mol, c, two_elec, compact=False)).reshape(
            (hcore.shape[0], hcore.shape[0], hcore.shape[0], hcore.shape[0]))
    }

    result['energy_nuc'] = mol.energy_nuc()
    result['energy_elec'] = mf.energy_elec()[0]

    return result

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
        second_q += FermionicOp(coeff, f'{k}^ {j}^ {i} {l}') + \
                    FermionicOp(coeff, f'{k}^ {j+n}^ {i+n} {l}') + \
                    FermionicOp(coeff, f'{k+n}^ {j}^ {i} {l+n}') + \
                    FermionicOp(coeff, f'{k+n}^ {j+n}^ {i+n} {l+n}')
    return second_q.jordan_wigner
