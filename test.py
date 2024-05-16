from qc_practice.fermionic_sort import FermionicSort
from qc_practice.jordan_wigner import JordanWigner

# operator = input("Input: ")
OPERATOR = '-1.00 0^ 0 -0.5 1^ 1 -0.3 2^ 2 +0.3 0^ 0^ 0 0 +0.3 1^ 1^ 1 1 +0.2 0^ 1^ 0 1 +0.2 1^ 2^ 2 1'
fermi_sort = FermionicSort(OPERATOR)
print("\nFermionic Sort:", fermi_sort)

jw = JordanWigner(fermi_sort)

print("Jordan Wigner Sort:", jw.pauli_strings)
