from itertools import product
from jqc.mapper.fermion import FermionicOp

def dipole_moment(profile):
    'Generate the dipole moment in the second quantization form.'
    no = profile.num_orb
    mo_dip = profile.qm['dipole']
    second_q = [FermionicOp() for _ in range(3)]
    for xyz in range(3):
        for i, j in product(range(no), repeat=2):
            dipole_term = mo_dip[xyz, i, j]
            second_q[xyz] += FermionicOp(dipole_term, f'{j}^ {i}') + \
                             FermionicOp(dipole_term, f'{i + no}^ {j + no}')
    return [sq.jordan_wigner for sq in second_q]
