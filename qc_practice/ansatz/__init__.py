# __init__.py

from typing import Protocol
from scipy.optimize import OptimizeResult
from .uccsd import UCCSD, UpCCSD, UCCGSD, UpCCGSD, kUpCCGSD
from .hea import HEA
from .sp import SP, RSP, OSP

class Ansatz(Protocol):
    'Protocol for the ansatz class.'
    def generate_coeff(self, profile, coeff: float = 1e-5) -> list[float]:
        'Generate ansatz coefficients.'
        ...
    def ansatz(self, qc, profile, coeff) -> None:
        'The ansatz for the calculation.'
        ...
    @staticmethod
    def call_optimizer(func, coeff, method) -> OptimizeResult:
        'Optimize the coefficients'
        ...
