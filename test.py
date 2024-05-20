from qc_practice.fermionic_sort import FermionicSort
from qc_practice.jordan_wigner import JordanWigner

# operator = input("Input: ")
OPERATOR = '''
+ 0.3378550774017582 0^ 0^ 0 0
+ 0.3322908651276484 0^ 1^ 1 0
+ 0.3378550774017582 0^ 2^ 2 0
+ 0.0904655998921157 0^ 3^ 3 0'''
# + 0.3378550774017582 2^ 0^ 0 2
# + 0.3378550774017582 2^ 2^ 2 2
# + 0.3322908651276484 0^ 3^ 2 1
# + 0.3322908651276484 2^ 1^ 0 3
# + 0.3322908651276484 2^ 3^ 3 2
# + 0.0904655998921157 0^ 0^ 1 1
# + 0.0904655998921157 0^ 2^ 3 1
# + 0.0904655998921157 2^ 0^ 1 3
# + 0.0904655998921157 2^ 2^ 3 3
# + 0.0904655998921157 0^ 1^ 0 1
# + 0.0904655998921157 2^ 1^ 1 2
# + 0.0904655998921157 2^ 3^ 2 3
# + 0.0904655998921157 1^ 0^ 1 0
# + 0.0904655998921157 1^ 2^ 2 1
# + 0.0904655998921157 3^ 0^ 0 3
# + 0.0904655998921157 3^ 2^ 3 2
# + 0.0904655998921157 1^ 1^ 0 0
# + 0.0904655998921157 1^ 3^ 2 0
# + 0.0904655998921157 3^ 1^ 0 2
# + 0.0904655998921157 3^ 3^ 2 2
# + 0.3322908651276482 1^ 0^ 0 1
# + 0.3322908651276482 1^ 2^ 3 0
# + 0.3322908651276482 3^ 0^ 1 2
# + 0.3322908651276482 3^ 2^ 2 3
# + 0.3492868613660088 1^ 1^ 1 1
# + 0.3492868613660088 1^ 3^ 3 1
# + 0.3492868613660088 3^ 1^ 1 3
# + 0.3492868613660088 3^ 3^ 3 3'''
fermi_sort = FermionicSort(OPERATOR)
print("\nFermionic Sort:", fermi_sort)

jw = JordanWigner(fermi_sort)

print("Jordan Wigner Sort:")
print(jw)

#%%

from qc_practice.pauli_string import PauliString
from qc_practice.pauli_operator import PauliOperator

string1 = PauliString([PauliOperator('Z'), PauliOperator('X')])
print(string1)