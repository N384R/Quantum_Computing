import numpy as np

class Pauli:
    def __init__(self, pauli):
        if pauli not in ('X', 'Y', 'Z', 
                         'iX', 'iY', 'iZ', 
                         '-X', '-Y', '-Z', 
                         '-iX', '-iY', '-iZ',
                         'I', '-I', 'iI', '-iI'):
            print("Error: Invalid Pauli Operator")
            exit()
        self.symbol = pauli

    def __mul__(self, other):
        if not isinstance(other, Pauli):
            print("Error: Invalid Pauli Operator")
            exit()
        if self.symbol == other.symbol and 'i' not in self.symbol:
            return Pauli('I')

        pauli, sign = self.calculation(self.symbol, other.symbol)
        return Pauli(sign + pauli)

    @property
    def matrix(self):
        'The matrix representation of the Pauli operator.'
        if self.symbol == 'X':
            return np.array([[0, 1], [1, 0]])
        elif self.symbol == 'Y':
            return np.array([[0, -1j], [1j, 0]])
        elif self.symbol == 'Z':
            return np.array([[1, 0], [0, -1]])
        elif self.symbol == 'I':
            return np.array([[1, 0], [0, 1]])

    @staticmethod
    def calculation(pauli1, pauli2):
        rules = {
            ('X', 'Y'): ('Z', 'i'),
            ('Y', 'Z'): ('X', 'i'),
            ('Z', 'X'): ('Y', 'i'),
            ('Y', 'X'): ('Z', '-i'),
            ('Z', 'Y'): ('X', '-i'),
            ('X', 'Z'): ('Y', '-i'),
            ('X', 'I'): ('X', ''),
            ('Y', 'I'): ('Y', ''),
            ('Z', 'I'): ('Z', ''),
        }

        sign1, op1 = Pauli.extract(pauli1)
        sign2, op2 = Pauli.extract(pauli2)

        if (op1, op2) in rules:
            _op, sign3 = rules[(op1, op2)]
        elif (op2, op1) in rules:
            _op, sign3 = rules[(op2, op1)]
        else:
            _op = 'I'
            sign3 = ''

        _sign = Pauli.calc_sign(sign1, sign2, sign3)
        return _op, _sign

    @staticmethod
    def extract(pauli):
        sign = ''
        if pauli.startswith('-'):
            sign = '-'
            pauli = pauli[1:]
        if pauli.startswith('i'):
            sign += 'i'
            pauli = pauli[1:]

        return sign, pauli

    @staticmethod
    def calc_sign(sign1, sign2, sign3):
        count_m = sign1.count('-') + sign2.count('-') + sign3.count('-')
        count_i = sign1.count('i') + sign2.count('i') + sign3.count('i')

        sign = ''
        if count_m % 2 == 1:
            sign += '-'
        if count_i // 2 == 1:
            sign = '-' if sign == '' else ''
        if count_i % 2 == 1:
            sign += 'i'

        return sign

    def __repr__(self):
        return f'{self.symbol}'


if __name__ == "__main__":
    X = Pauli('X')
    Y = Pauli('Y')
    Z = Pauli('Z')
    iZ = Pauli('iZ')

    test = iZ * Z
    print(test)
