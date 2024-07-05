'''
This module contains the implementation of the SSVQE class. 
The class is initialized with the molecule and the ansatz to be used.

The class has the following methods:

__init__: Initializes the VQE class with the molecule and the ansatz to be used.
circuit : Builds the quantum circuit for the calculation.
batch : Performs the VQE calculation for a given set of coefficients.
run : Runs the VQE optimization.
'''

import datetime
from pyscf.gto import Mole
from qiskit import QuantumCircuit
from qc_practice.ansatz import Ansatz
from qc_practice.ansatz.uccsd import UpCCGSD
from qc_practice.measure.measure import measure
from qc_practice.measure.hamiltonian import hamiltonian
from qc_practice.simulator import Simulator
from qc_practice.simulator.qasm import QASM
from qc_practice.profile import Profile

class VQE:
    '''
    Variational Quantum Eigensolver (VQE).

    Examples:
    >>> from pyscf import gto
    >>> from qc_practice import VQE
    >>> from qc_practice.ansatz import UCCSD
    >>> from qc_practice.simulator import QASM
    >>> mol = gto.M(atom = 'H 0 0 0; H 0 0 0.7', basis = 'sto-3g')
    >>> vqe = VQE(mol, ansatz = UCCSD(), simulator = QASM())
    >>> result = vqe.run()
    '''

    def __init__(self, mol: Mole, ansatz: Ansatz = UpCCGSD(), simulator: Simulator = QASM()):
        self.ansatz = ansatz
        self.simulator = simulator
        self.profile: Profile = Profile(mol)
        self.hamiltonian = hamiltonian(self.profile)
        self._config = {'iteration': 0, 'optimizer': 'powell', 'parallel': False}

    @property
    def ansatz(self):
        'The ansatz to be used for the calculation.'
        return self.__ansatz

    @ansatz.setter
    def ansatz(self, ansatz):
        self.__ansatz = ansatz

    @property
    def simulator(self) -> Simulator:
        'The simulator to be used for the calculation.'
        return self.__simulator

    @simulator.setter
    def simulator(self, simulator):
        self.__simulator = simulator

    @property
    def optimizer(self) -> str:
        'The optimizer to be used for the calculation.'
        return self._config['optimizer']

    @optimizer.setter
    def optimizer(self, optimizer):
        self._config['optimizer'] = optimizer

    def iteration(self) -> int:
        'Increments the iteration count.'
        self._config['iteration'] += 1
        return self._config['iteration']

    def circuit(self, coeff):
        'Builds the quantum circuit for the calculation.'
        qc = QuantumCircuit(self.profile.num_orb*2, self.profile.num_orb*2)
        for i in range(self.profile.num_elec//2):
            qc.x(i)
            qc.x(i + self.profile.num_elec)
        self.ansatz.ansatz(qc, self.profile, coeff)
        return qc

    @staticmethod
    def batch_output(func):
        'Decorator for the batch output.'
        def wrapper(self, *args, **kwargs):
            energy = func(self, *args, **kwargs)
            print(f"Iteration: {self.iteration()}, Energy: {self.profile.energy_total():12.09f}",
                  end='\r', flush=True)
            return energy
        return wrapper

    @batch_output
    def batch(self, coeff):
        'Performs the calculation for a given set of coefficients.'
        qc = self.circuit(coeff)
        energy = measure(qc, self.hamiltonian, self.simulator)
        self.profile.energy_elec = energy
        self.profile.coeff = coeff
        self.profile.circuit = qc
        return energy

    @staticmethod
    def general_output(func):
        'Decorator for the normal output.'
        def wrapper(self, *args, **kwargs):
            print(f'\nStarting {self.__class__.__name__} Calculation\n')
            print(f'Ansatz: {self.ansatz.__class__.__name__}')
            print(f'Simulator: {self.simulator.__class__.__name__}')
            print(f'Optimizer: {self.optimizer}\n')
            start = datetime.datetime.now()
            result = func(self, *args, **kwargs)
            elapsed = str(datetime.datetime.now() - start)
            print(f'Elapsed time: {elapsed.split(".", maxsplit=1)[0]}')
            return result
        return wrapper

    @general_output
    def run(self):
        'Performs the VQE calculation.'
        return self._run()

    @staticmethod
    def normal_output(func):
        'Decorator for the normal output.'
        def wrapper(self, *args, **kwargs):
            print(f'State {self.profile.state}:')
            result = func(self, *args, **kwargs)
            print(f"Iteration: {self.iteration()}, Converged!!         ")
            print(f'Total Energy: {result.energy_total():12.09f}\n')
            return result
        return wrapper

    @normal_output
    def _run(self):
        coeff = self.ansatz.generate_coeff(self.profile)
        optimized = self.ansatz.call_optimizer(self.batch, coeff, self.optimizer)
        self.profile.energy_elec = optimized.fun
        self.profile.coeff = optimized.x
        self.profile.circuit = self.circuit(optimized.x)
        return self.profile
