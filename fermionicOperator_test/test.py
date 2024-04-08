from fermionic_sort import FermionicSort
from jordanwigner import JordanWigner

# operator = input("Input: ")
operator = '3.87 0 2 1^ 2^ + 1.5 0 1^'
fermi_sort = FermionicSort(operator)
print("Fermionic Sort:", fermi_sort)

jw = JordanWigner(fermi_sort)

print(jw)
