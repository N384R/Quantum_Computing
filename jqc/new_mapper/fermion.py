class Fermion:
    'Class for Fermionic operator.'
    def __init__(self, fermion):
        self.fermion = fermion

    @property
    def type(self):
        'Return the type of fermion.'
        return 'creation' if '^' in self.fermion else 'annihilation'

    @property
    def num(self):
        'Return the number of the fermion.'
        if self.type == 'creation':
            return int(self.fermion.replace('^', ''))
        return int(self.fermion)

    def __repr__(self):
        sub = str.maketrans("0123456789^", "₀₁₂₃₄₅₆₇₈₉†")
        notation = self.fermion.translate(sub)
        line = "a" + notation
        return line

class FermionString:
    'Class for string of Fermionic operators.'
    def __init__(self, string):
        self._operators = split_op(string)

    @property
    def operators(self):
        'Return the operators in the string.'
        return self._operators

    def __add__(self, other):
        if isinstance(other, FermionString):
            return FermionStrings(self, other)
        if isinstance(other, FermionStrings):
            return FermionStrings(self, *other)
        return NotImplemented

    def __repr__(self):
        return " ".join([str(op) for op in self.operators])

class FermionStrings:
    'multiple FermionString objects.'
    def __init__(self, *args):
        self.fermion_strings = args
    
    def 

def split_op(operators):
    'Split the operators in a string and return a list of Fermion objects.'
    operators = operators.split()
    return [Fermion(op) for op in operators]

if __name__ == "__main__":
    fermion_list = ['2', '1^', '3', '2^', '1']
    fermion_objects = [Fermion(fermion) for fermion in fermion_list]
    print(fermion_objects)
