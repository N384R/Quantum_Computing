'''
This module contains functions to calculate angular momentum operators.
'''

from jqc.mapper.fermion import FermionicOp

def s_z(profile):
    'Calculates the Sz operator.'
    no = profile.num_orb
    operator = FermionicOp()
    for i in range(no):
        operator += FermionicOp(0.5, f'{i}^ {i}') - \
                    FermionicOp(0.5, f'{i+no}^ {i+no}')
    return operator

def s_minus(profile):
    'Calculates the S- operator.'
    no = profile.num_orb
    operator = FermionicOp()
    for i in range(no):
        operator += FermionicOp(1.0, f'{i+no}^ {i}')
    return operator

def s_plus(profile):
    'Calculates the S+ operator.'
    no = profile.num_orb
    operator = FermionicOp()
    for i in range(no):
        operator += FermionicOp(1.0, f'{i}^ {i+no}')
    return operator

def total_spin(profile):
    'Calculates the total spin operator.'
    sz_op = s_z(profile)
    sp_op = s_plus(profile)
    sm_op = s_minus(profile)

    s2_op = ((sp_op * sm_op) + (sm_op * sp_op)) / 2 + (sz_op * sz_op)
    return s2_op.jordan_wigner
