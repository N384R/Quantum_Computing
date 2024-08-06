from copy import deepcopy
from jqc.new_mapper.fermion import Fermion

def split_op(operators):
    operators = operators.split()
    return [Fermion(op) for op in operators]
