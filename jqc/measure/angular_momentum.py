'''
This module contains functions to calculate angular momentum operators.
'''

from jqc.mapper.jordan_wigner import JordanWignerMapper

def s_z(profile):
    'Calculates the Sz operator.'
    no = profile.num_orb
    operator = ''
    for i in range(no):
        operator += f'+ 0.5 {i}^ {i} - 0.5 {i+no}^ {i+no}' + '\n'
    return JordanWignerMapper(operator)

def s_minus(profile):
    'Calculates the S- operator.'
    no = profile.num_orb
    operator = ''
    for i in range(no):
        operator += f'+ 1.0 {i+no}^ {i}' + '\n'
    return JordanWignerMapper(operator)

def s_plus(profile):
    'Calculates the S+ operator.'
    no = profile.num_orb
    operator = ''
    for i in range(no):
        operator += f'+ 1.0 {i}^ {i+no}' + '\n'
    return JordanWignerMapper(operator)
