import numpy as np

class PauliOperator:
    def __init__(self, pauli):
        self.symbol = pauli
        self.matrix = self.get_operator(pauli)
    
    def get_operator(self, pauli):
        if pauli == 'X':
            return np.array([[0, 1], [1, 0]])
        elif pauli == 'Y':
            return np.array([[0, -1j], [1j, 0]])
        elif pauli == 'Z':
            return np.array([[1, 0], [0, -1]])
        elif pauli == 'I':
            return np.array([[1, 0], [0, 1]])
        elif pauli == 'f^':
            return np.array([[0, 1], [0, 0]])
        elif pauli == 'f':
            return np.array([[0, 0], [1, 0]])
        else:
            raise ValueError('Invalid Pauli operator')
    
    def __add__(self, other):
        try:
            return self.matrix + other.matrix
        except:
            return self.matrix + other
        
    def __radd__(self, other):
        return other + self.matrix
    
    def __sub__(self, other):
        try:
            return self.matrix - other.matrix
        except:
            return self.matrix - other
    
    def __mul__(self, other):
        try:
            return np.dot(self.matrix, other.matrix)
        except:
            return self.matrix * other
    
    def __rmul__(self, other):
        return other * self.matrix

    def __repr__(self):
        return f'{self.symbol}'
    
if __name__ == "__main__":
    X = PauliOperator('X')
    Y = PauliOperator('Y')
    Z = PauliOperator('Z')
    
    print(X)
    print(Y)
    print(Z)