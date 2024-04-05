from fermionic_sort import FermionicSort
from jordanwigner import JordanWigner

operator = input("Input: ")
fermi_sort = FermionicSort(operator)
print("Fermionic Sort:", fermi_sort)
jw = JordanWigner(fermi_sort)
print(jw)