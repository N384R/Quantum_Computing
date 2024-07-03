# __init__.py

from typing import Protocol
from .simulator import QASM, StateVector

class Simulator(Protocol):
    'Protocol for running simulators.'
    def run_simulator(self, qc, p_string) -> float:
        'Run the simulator.'
        ...
    def swap_test(self, state1, state2, ansatz) -> float:
        'Swap test for measuring the overlap between two states.'
        ...
