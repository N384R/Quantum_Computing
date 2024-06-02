from fermion_operator import FermionOperator
from pauli_operator import PauliOperator

class JordanWigner:
    def __init__(self, fermion_operator, orbital_number):
        self.fermion_operator = fermion_operator
        self.orbital_number = orbital_number
        self.pauli_operator = self.JWtransformation(fermion_operator, orbital_number)
    
    def __repr__(self):
        return f'{self.pauli_operator}'

    def JWtransformation(self, fermion_operator, orbital_number):
        pauli_string1 = {}
        pauli_string2 = {}
        z_string = {}
        X = PauliOperator('X')
        Y = PauliOperator('Y')
        Z = PauliOperator('Z')

        if fermion_operator.symbol == 'creation':
            pauli_string1[orbital_number] = lambda: (X - 1j*Y)/2
        elif fermion_operator.symbol == 'annihilation':
            pauli_string1[orbital_number] = lambda: (X + 1j*Y)/2

        for i in range(orbital_number):
                z_string[i] = lambda: Z
        
        pauli_string = {**pauli_string1, **pauli_string2, **z_string}
        return pauli_string

if __name__ == "__main__":
    fermion = FermionOperator('creation')
    JW = JordanWigner(fermion, 6)
    for key in JW.pauli_operator:
        print(JW.pauli_operator[key]())