# __init__.py

from typing import Protocol
from qiskit import QuantumCircuit
from qc_practice.profile import Profile
from .qasm import QASM
from .state_vector import StateVector

class Simulator(Protocol):
    'Protocol for running simulators.'
    def measure(self, qc: QuantumCircuit, operator, parallel: bool) -> float:
        'Measure the expectation value of a Hamiltonian'
        ...
    def get_overlap(self, state1: Profile, state2: Profile) -> float:
        'Get the square of the overlap between two states.'
        ...
