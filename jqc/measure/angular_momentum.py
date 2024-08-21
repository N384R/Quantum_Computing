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
    return operator.jordan_wigner

def s_minus(profile):
    'Calculates the S- operator.'
    no = profile.num_orb
    operator = FermionicOp()
    for i in range(no):
        operator += FermionicOp(1.0, f'{i+no}^ {i}')
    return operator.jordan_wigner

def s_plus(profile):
    'Calculates the S+ operator.'
    no = profile.num_orb
    operator = FermionicOp()
    for i in range(no):
        operator += FermionicOp(1.0, f'{i}^ {i+no}')
    return operator.jordan_wigner
