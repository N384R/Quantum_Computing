# __init__.py

from typing import Protocol
from qiskit import QuantumCircuit
from jqc.vqe.profile import Profile
from .qasm import QASM
from .state_vector import StateVector


class Simulator(Protocol):
    'Protocol for running simulators.'

    @property
    def parallel(self) -> bool:
        'The parallel flag for the calculation.'
        ...

    @parallel.setter
    def parallel(self, value: bool):
        ...

    def measure(self, qc1, operator, qc2=None) -> float:
        'Measure the expectation value of a Hamiltonian'
        ...

    def get_overlap(self, state1: Profile, state2: Profile) -> float:
        'Get the square of the overlap between two states.'
        ...
