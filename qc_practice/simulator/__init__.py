# __init__.py

from typing import Protocol
from qc_practice.profile import Profile
from .qasm import QASM
from .state_vector import StateVector

class Simulator(Protocol):
    'Protocol for running simulators.'
    def run_simulator(self, qc, p_string) -> float:
        'Run the simulator.'
        ...
    def get_overlap(self, state1: Profile, state2: Profile) -> float:
        'Get the square of the overlap between two states.'
        ...
