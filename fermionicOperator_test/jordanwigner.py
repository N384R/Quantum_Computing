from fermion import Fermion
from pauli_operator import PauliOperator
from pauli_string import PauliString, PauliStrings
from fermionic_sort import FermionicSort

class _JordanWigner:
    def __init__(self, fermionstring, maximum):
        self.fermionstring = [Fermion(fermion) for fermion in fermionstring]
        self.maximum = maximum
        self.paulistrings = self.pauli_convert(self.fermionstring)
        
    def __iter__(self):
        return iter(self.paulistrings)

    def _pauli_convert(self, fermion):
        paulistring1 = PauliString()
        paulistring2 = PauliString()
        z_string = PauliString()
        i_string = PauliString()

        num = abs(fermion.num)
        for i in range(num):
            z_string[i] = PauliOperator('Z')
        
        if num < self.maximum:
            for i in range(num+1, self.maximum+1):
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

    def pauli_convert(self, fermionstring):
        result = None
        for fermion in fermionstring:
            pauli = self._pauli_convert(fermion)
            if result is None:
                result = pauli
            else:
                result *= pauli
        return result
    
    def __repr__(self):
        return f'{self.paulistrings}'

class JordanWigner():
    def __init__(self, fermionstrings):
        self.fermionstrings = fermionstrings
        self.maximum = self.max_num(fermionstrings)
        self.paulistrings = self.call_jordan_wigner(self.fermionstrings)

    def call_jordan_wigner(self, fermionstrings):
        result = []
        for fermionstring in fermionstrings:
            if fermionstring in ('-', '+'):
                sign = fermionstring
                continue
            
            coeff = float(fermionstring[0])/(2**(self.maximum+1))
            pauli = _JordanWigner(fermionstring[1:], self.maximum)
            for p in pauli:
                p[0] *= PauliOperator('-I') if sign == '-' else PauliOperator('I')
                result.append([coeff, p])
        return result
    
    def max_num(self, fermionstrings):
        maximum = 0
        for fermionstring in fermionstrings:
            if fermionstring in ('-', '+'):
                continue
            num = max([Fermion(fermion).num for fermion in fermionstring[1:]])
            maximum = num if num > maximum else maximum
        return maximum

    def __repr__(self):
        line = ''
        for i in range(len(self.paulistrings)):
            symbol = ''.join([f'{val}' for val in self.paulistrings[i][1].values()])
            count_m, count_i = symbol.count('-'), symbol.count('i')
            line += '-' if count_m == 1 else '+'
            line += f' ({self.paulistrings[i][0]})'
            line += 'i' if count_i == 1 else ''
            line += ' ' + symbol.replace('-', '').replace('i', '') + ' '
        return line


if __name__ == '__main__':
    jw = JordanWigner(['-', ['2.54', '1^', '3^', '2', '0']])
    print(jw)