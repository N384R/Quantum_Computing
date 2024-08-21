'''
This module contains Pauli class and PauliOp class.
'''

from functools import reduce
from itertools import chain
import operator
import numpy as np

class Pauli:
    'A class to represent a Pauli operator.'

    valid_paulis = ('X', 'iX', '-X', '-iX',
                    'Y', 'iY', '-Y', '-iY',
                    'Z', 'iZ', '-Z', '-iZ',
                    'I', 'iI', '-I', '-iI')

    rules = {('X', 'I'):   'X', ('Y', 'I'):   'Y', ('Z', 'I'):   'Z',
             ('I', 'X'):   'X', ('I', 'Y'):   'Y', ('I', 'Z'):   'Z',
             ('X', 'Y'):  'iZ', ('Y', 'Z'):  'iX', ('Z', 'X'):  'iY',
             ('Y', 'X'): '-iZ', ('Z', 'Y'): '-iX', ('X', 'Z'): '-iY'}

    matrices = {'X': np.array([[0,   1], [1,  0]]),
                'Y': np.array([[0, -1j], [1j, 0]]),
                'Z': np.array([[1,   0], [0, -1]]),
                'I': np.array([[1,   0], [0,  1]])}

    def __init__(self, pauli):
        if pauli not in Pauli.valid_paulis:
            raise ValueError("Invalid Pauli Operator")
        self._symbol = pauli

    @property
    def symbol(self):
        'Return the symbol of the Pauli operator.'
        return self._symbol

    @property
    def attribute(self):
        'Return the attribute of the Pauli operator.'
        sign = -1 if '-' in self.symbol else 1
        comp = 1j if 'i' in self.symbol else 1
        return (sign, comp, self.symbol[-1])

    @property
    def matrix(self):
        'Return the matrix representation of the Pauli operator.'
        sign, comp, op = self.attribute
        matrix = Pauli.matrices[op]
        return sign * comp * matrix

    def __mul__(self, other):
        if not isinstance(other, Pauli):
            raise ValueError("Invalid Pauli Operator")
        return self.multiply(self, other)

    def __hash__(self):
        return hash(self.symbol)

    def __eq__(self, other):
        if isinstance(other, Pauli):
            return self.symbol == other.symbol
        return False

    @staticmethod
    def multiply(p1, p2):
        'Return the product of two Pauli operators.'
        *sgn1, op1 = p1.attribute
        *sgn2, op2 = p2.attribute

        sym = 'I' if op1 == op2 else Pauli.rules[(op1, op2)]
        *sgn3, op3 = Pauli(sym).attribute

        product = reduce(operator.mul, sgn1 + sgn2 + sgn3)
        sign = '-' if product.real < 0 or product.imag < 0 else ''
        comp = 'i' if product.imag != 0 else ''
        return Pauli(f'{sign}{comp}{op3}')

    def __repr__(self):
        return self.symbol

class PauliOp:
    'A class to represent a string of Pauli operators.'

    def __init__(self, *args):
        if not args:
            self._objects = {}
        elif isinstance(args[0], dict):
            self._objects = args[0]
        elif isinstance(args[1], str):
            self._objects = get_op(args[1], args[0])
        else:
            raise TypeError('Invalid input.')

    @property
    def objects(self):
        'Return the operators in the string.'
        return self._objects

    def __add__(self, other):
        result = self.objects.copy()
        if not isinstance(other, PauliOp):
            return NotImplemented
        for k, v in other.objects.items():
            result[k] = result.get(k, 0) + v
            if abs(result[k]) < 1e-15:
                del result[k]
        return PauliOp(result)

    def __sub__(self, other):
        result = self.objects.copy()
        if not isinstance(other, PauliOp):
            return NotImplemented
        for k, v in other.objects.items():
            result[k] = result.get(k, 0) - v
            if abs(result[k]) < 1e-15:
                del result[k]
        return PauliOp(result)

    def __mul__(self, other):
        if not isinstance(other, PauliOp):
            return NotImplemented
        result = {}
        for k, v in distribute_ops(self.objects, other.objects):
            result[k] = result.get(k, 0) + v
            if abs(result[k].real) < 1e-15:
                result[k] = result[k].imag * 1j
            if abs(result[k].imag) < 1e-15:
                result[k] = result[k].real
        return PauliOp(result)

    def __truediv__(self, other):
        if not isinstance(other, (int, float)):
            return NotImplemented
        result = {}
        for k, v in self.objects.items():
            result[k] = v / other
        return PauliOp(result)

    def items(self):
        'Return the items in the string.'
        return self.objects.items()

    def __repr__(self):
        def sgn(v):
            return '-' if v < 0 else '+'
        result = ''
        for pauli, val in self.objects.items():
            ops  = ''.join(str(op) for op in pauli)
            val_str = f'{sgn(val.real)} {abs(val.real):.06f} ' + \
                      f'{sgn(val.imag)} {abs(val.imag):.06f} i'
            result += f'+ ({val_str}) {ops}\n'
        return result

def get_op(obj, val) -> dict:
    'Return the operators in a string.'
    obj = tuple(Pauli(v) for v in obj.split())
    return {obj: val}

def arrange_ops(obj, val):
    'Return the reduced operators in a string.'
    *sgns, ops = zip(*[op.attribute for op in obj])
    sgns_flat  = list(chain(*sgns))
    product    = reduce(operator.mul, sgns_flat)
    ops_arr    = tuple(Pauli(op) for op in ops)
    val_arr    = val * product
    yield ops_arr, val_arr

def distribute_ops(obj1, obj2):
    'Return the product of two Pauli strings.'
    for key1, val in obj1.items():
        for key2 in obj2:
            op = pauli_mult(key1, key2)
            yield from arrange_ops(op, val)

def pauli_mult(obj1, obj2):
    'Return the product of two Pauli operators.'
    return tuple(op1 * op2 for op1, op2 in zip(obj1, obj2))

if __name__ == "__main__":
    pauli1 = Pauli('iX')
    pauli2 = Pauli('-iZ')
    pauli3 = pauli1 * pauli2
    print(pauli3)

    ps1 = PauliOp(1.00, 'Z Z')
    ps2 = PauliOp(1.00, 'X X')
    ps3 = ps1 * ps2
    print(ps3)
