from copy import deepcopy
from jqc.new_mapper.fermion import Fermion

def split_op(operators):
    operators = operators.split()
    string = (Fermion(op) for op in operators[1:])
    return {string: operators[0]}
