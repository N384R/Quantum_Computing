import numpy as np

class Dot:
    def __init__(self, A, B):
        self.result = np.dot(A.matrix, B.matrix)
    
    def __repr__(self):
        return f'{self.result}'
    
class Tensor:
    def __init__(self, A, B):
        self.result = np.kron(A.matrix, B.matrix)
    
    def __repr__(self):
        return f'{self.result}'