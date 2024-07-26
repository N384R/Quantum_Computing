from itertools import combinations
import numpy as np
from pyscf import ao2mo
from qc_practice.mapper.jordan_wigner import JordanWignerMapper
from qc_practice.ansatz.uccsd import sign_p

def dipole_moment(profile):
    'Generate the dipole moment in the second quantization form.'
    no = profile.num_orb
    coords = profile.qm['coords']
    second_q = {'x': '', 'y': '', 'z': ''}
    for i, j in combinations(range(len(coords)), 2):
        r_vector = coords[j] - coords[i]
        for idx, comp in enumerate('xyz'):
            r = r_vector[idx]
            second_q[comp] += f'{sign_p(r)} {abs(r)} {j}^ {i} '
            second_q[comp] += f'{sign_p(r)} {abs(r)} {j + no}^ {i + no} '

    return {'x': JordanWignerMapper(second_q['x']),
            'y': JordanWignerMapper(second_q['y']),
            'z': JordanWignerMapper(second_q['z'])}
