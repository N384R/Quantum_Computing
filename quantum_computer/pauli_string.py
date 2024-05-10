class PauliString:
    def __init__(self, pauli=None):
        self.paulistring = {}
        if pauli is not None:
            for key, value in enumerate(pauli):
                self.paulistring[key] = value
    
    def __getitem__(self, key):
        return self.paulistring[key]
    
    def __setitem__(self, key, value):
        self.paulistring[key] = value

    def __delitem__(self, key):
        del self.paulistring[key]

    def __iter__(self):
        return iter(self.paulistring)
    
    def __len__(self):
        return len(self.paulistring)
    
    def keys(self):
        return self.paulistring.keys()

    def values(self):
        return self.paulistring.values()
    
    def items(self):
        return self.paulistring.items()
    
    def __mul__(self, other):
        if not isinstance(other, PauliString):
            print("Error: Invalid Pauli Operator")
            exit()
        return self.string_calculation(self.paulistring, other.paulistring)

    def string_calculation(self, string1, string2):
        result = PauliString()
        for key in string1.keys():
            result[key] = string1[key] * string2[key]
        return result
    
    def __add__(self, other):
        if isinstance(other, PauliString):
            return PauliStrings(self, other)
        elif isinstance(other, PauliStrings):
            return PauliStrings(self, *other)

    def __repr__(self):
        line = ''.join([f'{val}' for val in self.paulistring.values()])
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
        self.paulistrings = sorted(args, key=lambda x: str(x).replace('i', ''))

    def __getitem__(self, key):
        return self.paulistrings[key]

    def __add__(self, other):
        if isinstance(other, PauliString):
            return PauliStrings(*self, other)
        elif isinstance(other, PauliStrings):
            return PauliStrings(*self, *other)

    def __mul__(self, other):
        result = PauliStrings()
        for string1 in self:
            for string2 in other:
                result += string1 * string2
        return result

    def __repr__(self):
        return ' '.join([f'{string}' for string in self.paulistrings])

if __name__ == '__main__':
    string1 = PauliString([PauliOperator('Z'), PauliOperator('iX'), 
                           PauliOperator('-Y'), PauliOperator('I')])
    print('1:', string1)

    string2 = PauliString([PauliOperator('Z'), PauliOperator('Z'),
                            PauliOperator('X'), PauliOperator('X')])
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