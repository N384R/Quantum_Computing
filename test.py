from qc_practice.fermionic_sort import FermionicSort
from qc_practice.jordanwigner import JordanWigner

# operator = input("Input: ")
OPERATOR = '3.87 0 2 1^ 2^ + 1.5 0 1^'
fermi_sort = FermionicSort(OPERATOR)
print("Fermionic Sort:", fermi_sort)

jw = JordanWigner(fermi_sort)

print("Jordan Wigner Sort:", jw)
