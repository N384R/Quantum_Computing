from fermion import Fermion
from pauli_operator import PauliOperator
from pauli_operator import PauliSort
from fermionic_sort import FermionicSort

class JordanWigner:
    def __init__(self, operator=None):
        self.operator = operator
        if operator is None:
            self.operator = input("Input: ")
        
        self.fermioperator = FermionicSort(self.operator)

    def PauliConvert(self, operator, maximum):
        fermion = Fermion(operator)
        pauli_string1 = {}
        pauli_string2 = {}
        z_string = {}
        i_string = {}
        for i in range(fermion.num):
            z_string[i] = PauliOperator('Z')
        
        pauli_string1[fermion.num] = PauliOperator('X')
        pauli_string2[fermion.num] = PauliOperator('-iY') if fermion.type == 'creation' else PauliOperator('iY')

        if fermion.num < maximum:
            for i in range(fermion.num+1, maximum+1):
                i_string[i] = PauliOperator('I')
        
        pauli_string1.update(z_string); pauli_string1.update(i_string)
        pauli_string2.update(z_string); pauli_string2.update(i_string)
        pauli_string1 = PauliSort(pauli_string1)
        pauli_string2 = PauliSort(pauli_string2)

        return pauli_string1, pauli_string2
        

if __name__ == '__main__':
    JW = JordanWigner("-1.00 1^ 2^ 3 4")
    print(JW.PauliConvert('1^', 3))