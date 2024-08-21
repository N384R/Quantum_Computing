import sys
import numpy as np
from .pauli_operator import Pauli

class PauliString:
    def __init__(self, pauli=None):
        self.pauli_string = {}
        if pauli is not None:
            for key, value in enumerate(pauli):
                self.pauli_string[key] = value
        self.symbol = self.get_symbol(self.pauli_string)

    def __getitem__(self, key):
        return self.pauli_string[key]

    def __setitem__(self, key, value):
        self.pauli_string[key] = value

    def __delitem__(self, key):
        del self.pauli_string[key]

    def __iter__(self):
        return iter(self.pauli_string)

    def __len__(self):
        return len(self.pauli_string)

    def keys(self):
        return self.pauli_string.keys()

    def values(self):
        return self.pauli_string.values()

    def items(self):
        return self.pauli_string.items()

    def get_symbol(self, pauli):
        return ''.join(f'{val}' for val in pauli.values())

    @property
    def matrix(self):
        'The matrix representation of the Pauli string.'
        matrix = np.array([[1]])
        for value in self.pauli_string.values():
            matrix = np.kron(value.matrix, matrix)
        return matrix

    @property
    def count_iden(self):
        'The number of identity operators in the Pauli string.'
        return sum(1 for val in self.pauli_string.values() if val.symbol == 'I')

    def __mul__(self, other):
        if not isinstance(other, PauliString):
            print("Error: Invalid Pauli Operator")
            sys.exit()
        return self._string_calculation(self.pauli_string, other.pauli_string)

    def _string_calculation(self, s1, s2):
        result = PauliString()
        for key in s1.keys():
            result[key] = s1[key] * s2[key]
        return result

    def __add__(self, other):
        if isinstance(other, PauliString):
            return PauliStrings(self, other)
        if isinstance(other, PauliStrings):
            return PauliStrings(self, *other)
        return NotImplemented

    def __eq__(self, other):
        if not isinstance(other, PauliString):
            return NotImplemented
        return self.symbol == other.symbol

    def __hash__(self):
        return hash(self.symbol)

    def __repr__(self):
        line = ''.join([f'{val}' for val in self.pauli_string.values()])
        count_m, count_i = line.count('-'), line.count('i')
        line = line.replace('-', '').replace('i', '')
        sign = ''
        if count_m % 2 == 1:
            sign += '- '
        if count_i // 2 == 1:
            sign = '- ' if sign == '' else ''
        if count_i % 2 == 1:
            sign += 'i'

        if '-' not in sign:
            sign = '+ ' + sign

        return sign + line

class PauliStrings:
    def __init__(self, *args):
        self.pauli_strings = args

    def __getitem__(self, key):
        return self.pauli_strings[key]

    def __add__(self, other):
        if isinstance(other, PauliString):
            return PauliStrings(*self, other)
        elif isinstance(other, PauliStrings):
            return PauliStrings(*self, *other)

    def __mul__(self, other):
        result = PauliStrings()
        for s1 in self:
            for s2 in other:
                result += s1 * s2
        return result

    def __repr__(self):
        return ' '.join([f'{string}' for string in self.pauli_strings])


if __name__ == '__main__':
    string1 = PauliString([Pauli('Z'), Pauli('iX'),
                           Pauli('-Y'), Pauli('I')])
    print('1:', string1)

    string2 = PauliString([Pauli('Z'), Pauli('Z'),
                            Pauli('X'), Pauli('X')])
    print('2:', string2)

    string3 = string1 * string2
    print('3:', string3)

    string4 = string1 + string2
    print('4:', string4)

    string5 = string1 + string3
    print('5:', string5)

    string6 = string4 + string5
    print('6:', string6)

    string7 = string4 * string5
    print('7:', string7.paulistrings)

    string8 = PauliString([Pauli('Z'), Pauli('Z'),
                            Pauli('X'), Pauli('X')])


    if string2 == string8:
        print(string2, string8, 'Equal')
    else:
        print(string2, string8, 'Not Equal')

    print(string2.pauli_string)
