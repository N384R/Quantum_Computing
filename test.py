from qc_practice.fermionic_sort import FermionicSort
from qc_practice.jw import JordanWigner

# operator = input("Input: ")
OPERATOR = '-1.120959456292115 0^ 0 -0.9593757718161824 0^ 1 -0.9593757718161824 1^ 0 -1.120959456292115 1^ 1'
fermi_sort = FermionicSort(OPERATOR)
print("\nFermionic Sort:", fermi_sort)

jw = JordanWigner(fermi_sort)

print("Jordan Wigner Sort:", jw)
