from copy import deepcopy
from fermion import Fermion

class Fermionic_operator:
    def __init__(self):
        self.operator = input("Input: ")
        self.operator = self.split_operator(self.operator)
        self.operator = self.number_sort(self.operator)
        self.operator = self.dirac_sort(self.operator)
        self.operator = self.number_sort(self.operator)
        self.operator = self.length_sort(self.operator)
        self.operator = self.compute_operator(self.operator)
        print(self.print_operator(self.operator))

    def split_operator(self, operator):
        operator = operator.split()
        _operator = []
        __operator = []
        if len(operator) == 0:
                print("Error: Invalid Operator")
                exit()

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
            if (operator[i] == '-') or (operator[i] == '+'):
                continue
            c, _operator = _number_sort(operator[i])
            if c and (_operator[0] == '-'):
                del _operator[0]
                operator[i] = _operator
                if operator[i-1] == '-':
                    operator.insert(i, '+')
                    del operator[i-1]
                elif operator[i-1] == '+':
                    operator.insert(i, '-')
                    del operator[i-1]
                else:
                    operator.insert(i, '-')
                return self.number_sort(operator)
        return operator
    
    def dirac_sort(self, operator):
        def _dirac_sort(operator, sign):
            _operator = []
            __operator = deepcopy(operator)
            fermion = [Fermion(fermion) for fermion in operator]
            for i in range(1, len(operator)-1):
                if (fermion[i].type == 'annihilation') and (fermion[i+1].type == 'creation'):
                    if operator[i] == operator[i+1][:-1]:
                        operator[i], operator[i+1] = operator[i+1], operator[i]
                        del __operator[i]; del __operator[i]
                        _operator.append(__operator)
                        if sign == '+':
                            _operator.append('-')
                        elif sign == '-':
                            _operator.append('+')
                        _operator.append(operator)
                        return True, True, _operator
                    else:
                        operator[i], operator[i+1] = operator[i+1], operator[i]
                        if operator[0] != '-':
                            operator.insert(0, '-')
                        elif operator[0] == '-':
                            del operator[0]
                        return True, False, operator
            return False, False, operator

        for i in range(len(operator)):
            if (operator[i] == '-') or (operator[i] == '+') or (len(operator[i]) == 1):
                continue
            a, b, _operator = _dirac_sort(operator[i], operator[i-1])
            if not a:
                continue
                
            if b:
                for j in reversed(range(len(_operator))):
                    operator.insert(i, _operator[j])
                del operator[i+len(_operator)-1]
                return self.dirac_sort(operator)
            else:
                operator[i] = _operator
                if operator[i][0] == '-':
                    del operator[i][0]
                    if operator[i-1] == '-':
                        operator[i-1] = '+'
                    elif operator[i-1] == '+':
                        operator[i-1] = '-'
                    else:
                        if operator[i-1] == '-':
                            operator[i-1] = '+'
                        else:
                            operator.insert(i, '-')
                return self.dirac_sort(operator)
        return operator

    def length_sort(self, operator):  
        length = []
        for i in range(len(operator)):
            if (operator[i] != '-') and (operator[i] != '+'):
                length.append([i, len(operator[i])])

        _length = deepcopy(length)
        _length.sort(key=lambda x: x[1])
        _operator = deepcopy(operator)
        for i in range(len(length)):
            operator[_length[i][0]-1] = _operator[length[i][0]-1]
            operator[_length[i][0]] = _operator[length[i][0]]
        return operator

    def compute_operator(self, operator):
        i = 0
        for i in range(len(operator)-2):
            try:
                if (operator[i] == '-') or (operator[i] == '+'):
                    continue
                if (operator[i][1:] == operator[i+2][1:]):
                    if operator[i-1] == '-':
                        coeff1 = - float(operator[i][0])
                    else:
                        coeff1 = float(operator[i][0])
                    if operator[i+1] == '-':
                        coeff2 = - float(operator[i+2][0])
                    else:
                        coeff2 = float(operator[i+2][0])
                    coeff = coeff1 + coeff2
                    operator[i][0] = coeff
                    del operator[i+1]; del operator[i+1]
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
        

fo = Fermionic_operator()
fo.operator