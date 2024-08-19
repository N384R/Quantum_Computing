'''
This module defines the Fermion class and FermionicOp class.
'''

from itertools import pairwise

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
        return -int(self.fermion)

    def __eq__(self, other) -> bool:
        if isinstance(other, Fermion):
            return self.fermion == other.fermion
        return False

    def __hash__(self) -> int:
        return hash(self.fermion)

    def __hash__(self):
        return hash(self.fermion)

    def __repr__(self):
        sub = str.maketrans("0123456789^", "₀₁₂₃₄₅₆₇₈₉†")
        notation = self.fermion.translate(sub)
        line = "a" + notation
        return line

class FermionicOp:
    'Class for string of Fermionic operators.'
    def __init__(self, *args):
        if isinstance(args[0], dict):
            self._objects = args[0]
        elif isinstance(args[0], float) and isinstance(args[1], str):
            self._objects = split_op(args[0], args[1])
        else:
            raise TypeError('Invalid input.')

    @property
    def objects(self):
        'Return the operators in the string.'
        return self._objects

    def __add__(self, other):
        result = self._objects.copy()
        if isinstance(other, FermionicOp):
            other = other._objects
        else:
            return NotImplemented
        for k, v in other.items():
            result[k] = result.get(k, 0) + v
        return FermionicOp(result)

    @staticmethod
    def split_op(operators):
        'Split the operators in the string.'
        operators = operators.split()
        string = tuple(Fermion(op) for op in operators[1:])
        return {string: float(operators[0])}

    def __repr__(self):
        result = ''
        for operator, value in self._objects.items():
            sign = '-' if value < 0 else '+'
            ops = ''.join(str(op) for op in operator)
            result += f'{sign} ({abs(value):.06f}) {ops}\n'
        return result

def split_op(val, obj) -> dict:
    'Split the operators in a string and return a list of Fermion objects.'
    obj = tuple(Fermion(v) for v in obj.split())
    return {op: v for v, op in fermionic_sort(val, obj)}

def fermionic_sort(val, obj):
    'Sort the operators in a string.'
    val, obj = num_sort(val, obj)
    for v, op in dirac_sort(val, obj):
        yield v, op

def num_sort(val, obj):
    'Sort the operators by the number of the fermion.'
    obj_list = list(obj)
    swapped = False
    for i, j in pairwise(range(len(obj_list))):
        if mode(obj_list[i], obj_list[j]) == 'swap':
            val, obj_list[i], obj_list[j] = -val, obj_list[j], obj_list[i]
            swapped = True
    if swapped:
        return num_sort(val, tuple(obj_list))
    return val, tuple(obj_list)

def mode(a, b):
    'Return the mode of sorting.'
    if a.type == b.type:
        if a.num > b.num:
            return 'swap'
    if b.type == 'creation':
        return 'dirac'
    return 'keep'

def dirac_sort(val, obj):
    'Sort the operators by the commutator relation.'
    obj_list = list(obj)
    for i, j in pairwise(range(len(obj_list))):
        if mode(obj_list[i], obj_list[j]) == 'dirac':
            if obj_list[i].num == -obj_list[j].num:
                reduced_op = obj_list[:i] + obj_list[j+1:]
                yield from dirac_sort(val, tuple(reduced_op))
            obj_list[i], obj_list[j] = obj_list[j], obj_list[i]
            yield from dirac_sort(-val, tuple(obj_list))
            return
    yield num_sort(val, tuple(obj_list))

if __name__ == "__main__":
    fo1 = FermionicOp(1.00, '1 2 2^ 1^')
    fo2 = FermionicOp(2.00, '1^ 2^ 2 1')
    fo3 = fo1 + fo2
    print(fo3)
