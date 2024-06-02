from spinor import Spinor
from fermion_operator import create, annihilate
from pauli_operator import PauliOperator
from calculator import Dot, Tensor

a = Spinor('1')
b = Spinor('0')

X = PauliOperator('X')
Y = PauliOperator('Y')
Z = PauliOperator('Z')

print('X * 1:\n',Dot(X, a))
print('\nY * 1:\n',Dot(Y, a))
print('\nZ * 1:\n',Dot(Z, a))

print('\nX * Y:\n', Dot(X, Y))
print('\nY * X:\n', Dot(Y, X))
print('\nX * Z:\n', Dot(X, Z))

print('\nX @ Y:\n', Tensor(X, Y))