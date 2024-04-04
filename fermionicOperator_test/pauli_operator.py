class PauliOperator:
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
        if not isinstance(other, PauliOperator):
            print("Error: Invalid Pauli Operator")
            exit()
        
        if self.symbol == other.symbol:
            return PauliOperator('I')

        pauli, sign = self.calculation(self.symbol, other.symbol)
        return PauliOperator(sign + pauli)
        
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

        sign1, op1 = PauliOperator.extract(pauli1)
        sign2, op2 = PauliOperator.extract(pauli2)

        if (op1, op2) in rules:
            _op, sign3 = rules[(op1, op2)]
        elif (op2, op1) in rules:
            _op, sign3 = rules[(op2, op1)]
        else:
            _op = 'I'
            sign3 = ''
    
        _sign = PauliOperator.calc_sign(sign1, sign2, sign3)
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
    
class PauliSort():
    def __init__(self, paulistring):
        self.paulistring = self.paulisort(paulistring)
    
    def paulisort(self, paulistring):
        sorted_string = sorted(paulistring.items(), key=lambda x: x[0])
        return dict(sorted_string)
    
    def __repr__(self):
        return ''.join([f'{val}' for val in self.paulistring.values()])
    
if __name__ == "__main__":
    X = PauliOperator('X')
    Y = PauliOperator('Y')
    Z = PauliOperator('Z')
    iZ = PauliOperator('iZ')

    test = iZ * Z
    print(test)

    # string = {'0': X, '2': X, '1': Z}
    # string = PauliSort(string)
    # print(string)