from copy import deepcopy
from .fermion import Fermion

class FermionicSort:
    def __init__(self, operator=None):
        self.operator = operator
        if operator is None:
            self.operator = input("Input: ")
        split_operator = self.split_operator(self.operator)
        number_sort = self.number_sort(split_operator)
        dirac_sort = self.dirac_sort(number_sort)
        number_sort = self.number_sort(dirac_sort)
        length_sort = self.length_sort(number_sort)
        compute_operator = self.compute_operator(length_sort)
        self.operator = compute_operator

    def __iter__(self):
        return iter(self.operator)

    def split_operator(self, operator):
        operator = operator.split()
        if not operator:
            print("Error: Invalid Operator")
            exit()

        for i, op in enumerate(operator):
            if any(sign in op for sign in ('-', '+')) and len(op) > 1:
                sign, operator[i] = operator[i][0], operator[i][1:]
                operator.insert(i, sign)

        result = [] 
        _operator = []
        for i, op in enumerate(operator):
            if op in ('-', '+'):
                if i != 0:
                    result.append(_operator)
                result.append(operator[i])
                _operator = []
                continue
            _operator.append(operator[i])
        result.append(_operator)

        return result

    def number_sort(self, operator):
        def _number_sort(operator):
            fermion = [Fermion(fermion) for fermion in operator]
            for i in range(1, len(operator)-1):
                if fermion[i].type == fermion[i+1].type:
                    if fermion[i].type == 'creation' and fermion[i].num > fermion[i+1].num:
                        operator[i], operator[i+1] = operator[i+1], operator[i]
                        operator.insert(0, '-')
                        return True, operator
                    elif fermion[i].type == 'annihilation' and fermion[i+1].num > fermion[i].num:
                        operator[i], operator[i+1] = operator[i+1], operator[i]
                        operator.insert(0, '-')
                        return True, operator
            return False, operator

        for i, op in enumerate(operator):
            if op in ('-', '+'):
                continue
            c, _operator = _number_sort(operator[i])
            if c and (_operator[0] == '-'):
                del _operator[0]
                operator[i] = _operator
                if operator[i-1] in ('-', '+'):
                    operator[i-1] = '+' if operator[i-1] == '-' else '-'
                else:
                    operator.insert(i, '-')
                return self.number_sort(operator)
        return operator

    def dirac_sort(self, operator):
        def _dirac_sort(operator, sign):
            _operator = deepcopy(operator)
            fermion = [Fermion(fermion) for fermion in operator]
            for i in range(1, len(operator)-1):
                if (fermion[i].type == 'annihilation') and (fermion[i+1].type == 'creation'):
                    if operator[i] == operator[i+1][:-1]:
                        operator[i], operator[i+1] = operator[i+1], operator[i]
                        del _operator[i:i+2]
                        return True, True, [_operator] + ['+' if sign == '-' else '-'] + [operator]
                    else:
                        operator[i], operator[i+1] = operator[i+1], operator[i]
                        if operator[0] != '-':
                            operator.insert(0, '-')
                        else:
                            operator.pop(0)
                        return True, False, operator
            return False, False, operator

        for i, op in enumerate(operator):
            if op in ('-', '+') or not op:
                continue

            a, b, _operator = _dirac_sort(operator[i], operator[i-1])
            if not a:
                continue
            if b:
                operator[i:i+1] = _operator
                return self.dirac_sort(operator)
            operator[i] = _operator
            if operator[i][0] == '-':
                del operator[i][0]
                if operator[i-1] in ('-', '+'):
                    operator[i-1] = '+' if operator[i-1] == '-' else '-'
                else:
                    operator.insert(i, '-')
            return self.dirac_sort(operator)
        return operator

    def length_sort(self, operator):  
        length = [
            [i, len(operator[i])]
            for i in range(len(operator)) 
            if operator[i] not in ('-', '+')
        ]
        _idx = deepcopy(length)
        _operator = deepcopy(operator)
        _idx.sort(key=lambda x: x[1])
        for i, idx in enumerate(length):
            operator[idx[0]], operator[idx[0]-1] = _operator[_idx[i][0]], _operator[_idx[i][0]-1]
        return operator

    def compute_operator(self, operator):
        _terms = {}
        _operator = []
        sign = '+'

        for item in operator:
            if isinstance(item, list):
                coeff, *terms = item
                terms_tuple = tuple(terms)
                coeff = float(coeff) if sign == '+' else -float(coeff)
                if terms_tuple in _terms:
                    _terms[terms_tuple] += coeff
                else:
                    _terms[terms_tuple] = coeff
            else:
                sign = item

        for terms, coeff in _terms.items():
            if coeff != 0:
                _operator.append('+' if coeff > 0 else '-')
                _operator.append([f'{abs(coeff)}'] + list(terms))
        return _operator

    def print_operator(self, operator):
        sub = str.maketrans("0123456789^", "₀₁₂₃₄₅₆₇₈₉†")
        line = ''
        for i, op in enumerate(operator):
            if op in ('-', '+'):
                line += ' ' + op
                continue

            _operator = deepcopy(operator[i])
            line += f' {_operator[0]:.16f}'
            for j in range(1, len(_operator)):
                notation = _operator[j].translate(sub)
                line += " a" + notation
        return line[1:]

    def __repr__(self):
        return f'{self.print_operator(self.operator)}'


if __name__ == "__main__":
    fo = FermionicSort('2.54 3^ 0^ 1 2')
    print('Output:', fo.operator)
    