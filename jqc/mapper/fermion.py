'''
This module defines the Fermion class and FermionicOp class.
'''

from itertools import pairwise
from .pauli import PauliOp


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
        if self.fermion.isdigit():
            return -int(self.fermion)
        return ValueError('Invalid fermion.')

    def __eq__(self, other) -> bool:
        if isinstance(other, Fermion):
            return self.fermion == other.fermion
        return False

    def __hash__(self) -> int:
        return hash(self.fermion)

    def __repr__(self):
        sub = str.maketrans("0123456789^", "₀₁₂₃₄₅₆₇₈₉†")
        notation = self.fermion.translate(sub)
        line = "a" + notation
        return line


class FermionicOp:
    'Class for string of Fermionic operators.'

    def __init__(self, *args):
        if not args:
            self._objects = {}
        elif isinstance(args[0], dict):
            self._objects = args[0]
        elif isinstance(args[0], (int, float)) and isinstance(args[1], str):
            self._objects = get_op(args[1], args[0])
        else:
            raise TypeError('Invalid input.')

    @property
    def objects(self):
        'Return the operators in the string.'
        return self._objects

    @property
    def jordan_wigner(self):
        'Convert a Fermionic operator to a Pauli operator.'
        return jordan_wigner(self.objects)

    def __add__(self, other):
        if not isinstance(other, FermionicOp):
            return NotImplemented
        result = self.objects.copy()
        for k, v in other.objects.items():
            result[k] = result.get(k, 0) + v
            if abs(result[k]) < 1e-12:
                del result[k]
        return FermionicOp(result)

    def __sub__(self, other):
        if not isinstance(other, FermionicOp):
            return NotImplemented
        result = self.objects.copy()
        for k, v in other.objects.items():
            result[k] = result.get(k, 0) - v
            if abs(result[k]) < 1e-12:
                del result[k]
        return FermionicOp(result)

    def __mul__(self, other):
        if isinstance(other, (int, float)):
            result = {}
            for k, v in self.objects.items():
                result[k] = v * other
            return FermionicOp(result)
        if isinstance(other, FermionicOp):
            result = {}
            for k, v in distribute_ops(self.objects, other.objects):
                result[k] = result.get(k, 0) + v
                if abs(result[k]) < 1e-12:
                    result[k] = 0
            return FermionicOp(result)
        return NotImplemented

    def __truediv__(self, other):
        if not isinstance(other, (int, float)):
            return NotImplemented
        result = {}
        for k, v in self.objects.items():
            result[k] = v / other
        return FermionicOp(result)

    def __repr__(self):
        result = ''
        for operator, value in self.objects.items():
            sign = '-' if value < 0 else '+'
            ops = ''.join(str(op) for op in operator)
            result += f'{sign} ({abs(value):.06f}) {ops}\n'
        return result


def distribute_ops(obj1, obj2):
    'Distribute the operators in a string.'
    for op1, val1 in obj1.items():
        for op2, val2 in obj2.items():
            yield from fermionic_sort(tuple([*op1, *op2]), val1 * val2)


def get_op(obj, val) -> dict:
    'Return the operators in a string.'
    obj = tuple(Fermion(v) for v in obj.split())
    return dict(fermionic_sort(obj, val))


def fermionic_sort(obj, val):
    'Sort the operators in a string.'
    obj, val = num_sort(obj, val)
    for op, v in dirac_sort(obj, val):
        if abs(v) > 1e-12:
            yield op, v


def num_sort(obj, val):
    'Sort the operators by the number of the fermion.'
    obj_list = list(obj)
    swapped = False
    for i, j in pairwise(range(len(obj_list))):
        if mode(obj_list[i], obj_list[j]) == 'swap':
            val, obj_list[i], obj_list[j] = -val, obj_list[j], obj_list[i]
            swapped = True
    if swapped:
        return num_sort(tuple(obj_list), val)
    return tuple(obj_list), val


def mode(a, b):
    'Return the mode of sorting.'
    if a.type == b.type:
        if a.num > b.num:
            return 'swap'
    if (a.type, b.type) == ('annihilation', 'creation'):
        return 'dirac'
    return 'keep'


def dirac_sort(obj, val):
    'Sort the operators by the commutator relation.'
    obj_list = list(obj)
    for i, j in pairwise(range(len(obj_list))):
        if mode(obj_list[i], obj_list[j]) == 'dirac':
            if obj_list[i].num == -obj_list[j].num:
                reduced_op = obj_list[:i] + obj_list[j+1:]
                yield from dirac_sort(tuple(reduced_op), val)
            obj_list[i], obj_list[j] = obj_list[j], obj_list[i]
            yield from dirac_sort(tuple(obj_list), -val)
            return
    yield num_sort(tuple(obj_list), val)


def jordan_wigner(fermionic_op) -> PauliOp:
    'Convert a Fermionic operator to a Pauli operator.'
    result = PauliOp()
    max_len = max(
        (abs(fo.num) + 1 for op in fermionic_op for fo in op), default=0)
    for op, val in fermionic_op.items():
        pstr = PauliOp(val, ' '.join(['I'] * max_len))
        if not op:
            result += pstr
            continue
        for fo in op:
            z = ['Z'] * abs(fo.num)
            y = ['iY'] if fo.type == 'creation' else ['-iY']
            i = ['I'] * (max_len - abs(fo.num) - 1)
            fo2po = PauliOp(val, ' '.join(z + ['X'] + i)) + \
                PauliOp(val, ' '.join(z + y + i))
            pstr *= fo2po
        result += pstr / (2 ** (len(op)))
    for op, val in result.objects.items():
        if sum(1 for po in op if po.symbol == 'Z') % 2:
            result.objects[op] *= -1
    return result


if __name__ == "__main__":
    fo1 = FermionicOp(1.00, '1 2 2^ 1^')
    fo2 = FermionicOp(2.00, '1^ 2^ 2 1')
    fo3 = fo1 + fo2
    print(f'Text:\n{fo3}')
    print(f'Dict:\n{fo3.objects}')
