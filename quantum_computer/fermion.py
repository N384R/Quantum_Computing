class Fermion:
    def __init__(self, fermion):
        self.fermion = fermion
        self.type = self.determine_type()
        self.num = self.determine_number()
        
    def determine_type(self):
        return 'creation' if '^' in self.fermion else 'annihilation'
    
    def determine_number(self):
        try:
            if self.type == 'creation':
                return int(self.fermion.replace('^', ''))
            elif self.type == 'annihilation':
                return int(self.fermion)
        except:
            self.type == 'others'
            return self.fermion
    
    def __repr__(self):
        sub = str.maketrans("0123456789^", "₀₁₂₃₄₅₆₇₈₉†")
        notation = self.fermion.translate(sub)
        line = "a" + notation
        return line
    
    @staticmethod
    def sort_fermions(fermions):
        sorted_fermions = sorted(fermions, key=lambda x: (x.type != 'creation', int(x.num)))
        return sorted_fermions

if __name__ == "__main__":
    fermion_list = ['2', '1^', '3', '2^', '1']
    fermion_objects = [Fermion(fermion) for fermion in fermion_list]

    sorted_fermions = Fermion.sort_fermions(fermion_objects)
    print(sorted_fermions)
