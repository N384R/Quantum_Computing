import numpy as np
from spinor import Spinor

class FermionOperator:
    def __init__(self, symbol):
        self.symbol = symbol
        self.matrix = self.get_operator(symbol)

    def get_operator(self, symbol):
        if symbol == 'creation':
            return np.array([[0, 1], [0, 0]])
        elif symbol == 'annihilation':
            return np.array([[0, 0], [1, 0]])
        else:
            raise ValueError('Invalid Fermion operator')
    
    def __repr__(self):
        return f'{self.matrix}'

# class FermionOperator(Spinor):
#     def __init__(self, spinor):
#         super().__init__(spinor.symbol)
#         self.operator = None
#         self.matrix = spinor.vector
#         self.symbol = spinor.symbol

#     def operation(self):
#         if self.operator is not None:
#             self.matrix = self.operator.dot(self.matrix)
#             self.symbol = self.get_symbol(self.matrix)
#         return self.matrix
    
#     def __repr__(self):
#         return f'{self.operation()}'

# class create(FermionOperator):
#     def __init__(self, spinor):
#         super().__init__(spinor)
#         self.operator = np.array([[0, 1], [0, 0]])

# class annihilate(FermionOperator):
#     def __init__(self, spinor):
#         super().__init__(spinor)
#         self.operator = np.array([[0, 0], [1, 0]])

# if __name__ == "__main__":
#     up = Spinor('1')
#     down = Spinor('0')
    
#     c_up = create(up)
#     c_down = create(down)
#     print(c_up)
#     print(c_down)
