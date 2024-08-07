class Fermion:
    def __init__(self, fermion):
        self.fermion = fermion

    @property
    def type(self):
        'Return the type of fermion.'
        if '^' in self.fermion:
            return 'creation'
        return 'annihilation'

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

if __name__ == "__main__":
    fermion_list = ['2', '1^', '3', '2^', '1']
    fermion_objects = [Fermion(fermion) for fermion in fermion_list]
    print(fermion_objects)
