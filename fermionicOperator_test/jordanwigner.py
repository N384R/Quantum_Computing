from fermion import Fermion
from pauli_operator import PauliOperator
from pauli_string import PauliString, PauliStrings
from fermionic_sort import FermionicSort

class JordanWigner:
    def __init__(self, fermionstring=None):
        self.fermionstring = [Fermion(fermion) for fermion in fermionstring]
        self.paulistrings = self.PauliConvert(self.fermionstring)

    def jordan_wigner(self, fermionstrings):
        result = []
        for fermionstring in fermionstrings:
            

    def _PauliConvert(self, fermion, maximum):
        paulistring1 = PauliString()
        paulistring2 = PauliString()
        z_string = PauliString()
        i_string = PauliString()

        num = abs(fermion.num)
        for i in range(num):
            z_string[i] = PauliOperator('Z')
        
        if num < maximum:
            for i in range(num+1, maximum+1):
                i_string[i] = PauliOperator('I')

        for key, value in z_string.items():
            paulistring1[key] = value
            paulistring2[key] = value
        
        paulistring1[num] = PauliOperator('X')
        paulistring2[num] = PauliOperator('-iY') if fermion.type == 'creation' else PauliOperator('iY')

        for key, value in i_string.items():
            paulistring1[key] = value
            paulistring2[key] = value

        return PauliStrings(paulistring1, paulistring2)

    def PauliConvert(self, fermionstring):
        maximum = max(fermion.num for fermion in fermionstring)
        result = None
        for fermion in fermionstring:
            pauli = self._PauliConvert(fermion, maximum)
            print(pauli)
            if result is None:
                result = pauli
            else:
                result *= pauli
        return result
    
    def __repr__(self):
        return f'{self.paulistrings}'


if __name__ == '__main__':

    jw = JordanWigner(['1^', '1'])
    print(jw)