from itertools import combinations
import numpy as np
from pyscf import ao2mo
from jqc.mapper.fermion import FermionicOp

def dipole_moment(profile):
    'Generate the dipole moment in the second quantization form.'
    no = profile.num_orb
    coords = profile.qm['coords']
    second_q = {'x': FermionicOp(),
                'y': FermionicOp(),
                'z': FermionicOp()}
    for i, j in combinations(range(len(coords)), 2):
        r_vector = coords[j] - coords[i]
        for idx, comp in enumerate('xyz'):
            r = r_vector[idx]
            second_q[comp] += FermionicOp(r, f'{j}^ {i}')
            second_q[comp] += FermionicOp(r, f'{j + no}^ {i + no}')

    return {'x': second_q['x'].jordan_wigner,
            'y': second_q['y'].jordan_wigner,
            'z': second_q['z'].jordan_wigner}
