import numpy as np

class Spinor:
    def __init__(self, spinor):
        self.symbol = spinor
        self.matrix = self.get_vector(spinor)
    
    def get_vector(self, spinor):
        if spinor == '1':
            return np.array([[1], [0]])
        elif spinor == '0':
            return np.array([[0], [1]])
        elif spinor == 'empty':
            return np.array([[0], [0]])
        else:
            raise ValueError('Invalid spinor')
        
    def get_symbol(self, vector):
        if np.array_equal(vector, np.array([[1], [0]])):
            return '1'
        elif np.array_equal(vector, np.array([[0], [1]])):
            return '0' 
        elif np.array_equal(vector, np.array([[0], [0]])):
            return 'empty'
        else:
            return ValueError('Invalid vector')
    
    def __repr__(self):
        return f'{self.symbol}'

if __name__ == "__main__":
    up = Spinor('1')
    down = Spinor('0')
    
    print(up)
    print(down)
    print(Spinor('empty'))