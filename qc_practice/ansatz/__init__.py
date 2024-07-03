# __init__.py

from typing import Protocol
from .uccsd import UCCSD, UpCCSD, UCCGSD, UpCCGSD, kUpCCGSD
from .hea import HEA
from .sp import SP, RSP, OSP

class Ansatz(Protocol):
    'Protocol for the ansatz class.'
    def generate_coeff(self, profile, coeff: float = 1e-5) -> list:
        'Generate ansatz coefficients.'
        ...
    def ansatz(self, qc, profile, coeff) -> None:
        'The ansatz for the calculation.'
