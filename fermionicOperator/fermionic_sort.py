from copy import deepcopy
from fermion import Fermion

class Fermionic_operator:
    def __init__(self, operator=None):
        if operator is None:
            self.operator = input("Input: ")
        split_operator = self.split_operator(self.operator)
        number_sort = self.number_sort(split_operator)
        dirac_sort = self.dirac_sort(number_sort)
        number_sort = self.number_sort(dirac_sort)
        length_sort = self.length_sort(number_sort)
        compute_operator = self.compute_operator(length_sort)
        self.operator = compute_operator

    def split_operator(self, operator):
        operator = operator.split()
        if not operator:
                print("Error: Invalid Operator")
                exit()

        _operator = []
        __operator = []

        for i in range(len(operator)):
            if operator[i] == '-':
                __operator.append(_operator)
                __operator.append('-')
                _operator = []
                continue
            elif operator[i] == '+':
                __operator.append(_operator)
                __operator.append('+')
                _operator = []
                continue
            _operator.append(operator[i])
        __operator.append(_operator)

        if '-' in __operator[0][0]:
            __operator[0][0] = __operator[0][0][1:]
            __operator.insert(0, '-')
        return __operator
    
    def number_sort(self, operator):
        def _number_sort(operator):
            fermion = [Fermion(fermion) for fermion in operator]
            for i in range(1, len(operator)-1):
                if (fermion[i].type == fermion[i+1].type) and (fermion[i].num > fermion[i+1].num):
                    operator[i], operator[i+1] = operator[i+1], operator[i]
                    operator.insert(0, '-')
                    return True, operator
            return False, operator 
        
        for i in range(len(operator)):
            if (operator[i] == '-') or (operator[i] == '+'): continue
            c, _operator = _number_sort(operator[i])
            if c and (_operator[0] == '-'):
                del _operator[0]
                operator[i] = _operator
                operator[i-1] = '+' if operator[i-1] == '-' else '-'
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
                        return True, True, [_operator] + ['-' if sign == '+' else '+'] + [operator]
                    else:
                        operator[i], operator[i+1] = operator[i+1], operator[i]
                        operator.insert(0, '-') if operator[0] != '-' else operator.pop(0)
                        return True, False, operator
            return False, False, operator

        for i in range(len(operator)):
            if operator[i] in ('-', '+') or len(operator[i]) == 1: 
                continue

            a, b, _operator = _dirac_sort(operator[i], operator[i-1])
            if not a: 
                continue
            
            if b:
                operator[i:i+1] = _operator
                return self.dirac_sort(operator)
            else:
                operator[i] = _operator
                if operator[i][0] == '-':
                    del operator[i][0]
                    if operator[i-1] in ['-', '+']:
                        operator[i-1] = '-' if operator[i-1] == '+' else '+'
                    else:
                        operator.insert(i, '-')
                return self.dirac_sort(operator)
        return operator

    def length_sort(self, operator):  
        length = [[i, len(operator[i])] for i in range(len(operator)) if operator[i] not in ('-', '+')]
        _length = deepcopy(length)
        _length.sort(key=lambda x: x[1])
        for i in range(len(length)):
            operator[length[i][0]], operator[length[i][0]-1] = operator[_length[i][0]], operator[_length[i][0]-1]
        return operator

    def compute_operator(self, operator):
        for i in range(len(operator)-2):
            try:
                if (operator[i] == '-') or (operator[i] == '+'):
                    continue
                if (operator[i][1:] == operator[i+2][1:]):
                    coeff1 = -float(operator[i][0]) if operator[i-1] == '-' else float(operator[i][0])
                    coeff2 = -float(operator[i+2][0]) if operator[i+1] == '-' else float(operator[i+2][0])
                    operator[i][0] = coeff1 + coeff2
                    del operator[i+2]
            except:
                break
        return operator

    def print_operator(self, operator):
        sub = str.maketrans("0123456789^", "₀₁₂₃₄₅₆₇₈₉†")
        
        line = ''
        for i in range(len(operator)):
            if (operator[i] == '-') or (operator[i] == '+'):
                line += ' ' + operator[i]
                continue

            _operator = deepcopy(operator[i])
            line += ' ' + str(_operator[0])
            for j in range(1, len(_operator)):
                notation = _operator[j].translate(sub)
                line += " a" + notation

        return line[1:]
        

if __name__ == "__main__":
    fo = Fermionic_operator()
    print('Output:', fo.print_operator(fo.operator))

    