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
from jqc.ansatz import Ansatz
from jqc.ansatz import UCCSD
from jqc.measure.hamiltonian import hamiltonian
from jqc.measure.angular_momentum import s_z
from jqc.simulator import Simulator
from jqc.simulator import StateVector
from jqc.vqe.profile import Profile

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

    def __init__(self, mol: Mole, ansatz: Ansatz = UCCSD(), simulator: Simulator = StateVector()):
        self.ansatz = ansatz
        self.simulator = simulator
        self.profile: Profile = Profile(mol)
        self.hamiltonian = hamiltonian(self.profile)
        self._config = {'iteration': 0, 'optimizer': 'powell', 'parallel': False, 'verbose': True}

    @property
    def ansatz(self):
        'The ansatz to be used for the calculation.'
        return self.__ansatz

    @ansatz.setter
    def ansatz(self, ansatz: Ansatz):
        self.__ansatz = ansatz

    @property
    def simulator(self):
        'The simulator to be used for the calculation.'
        return self.__simulator

    @simulator.setter
    def simulator(self, simulator: Simulator):
        self.__simulator = simulator

    @property
    def optimizer(self) -> str:
        'The optimizer to be used for the calculation.'
        return self._config['optimizer']

    @optimizer.setter
    def optimizer(self, optimizer: str):
        self._config['optimizer'] = optimizer

    @property
    def parallel(self) -> bool:
        'The parallel flag for the calculation.'
        return self._config['parallel']

    @parallel.setter
    def parallel(self, parallel: bool):
        self._config['parallel'] = parallel

    @property
    def verbose(self) -> bool:
        'The verbose flag for the calculation.'
        return self._config['verbose']

    @verbose.setter
    def verbose(self, verbose: bool):
        self._config['verbose'] = verbose

    def iteration(self) -> int:
        'Increments the iteration count.'
        self._config['iteration'] += 1
        return self._config['iteration']

    def circuit(self, coeff):
        'Builds the quantum circuit for the calculation.'
        qc = QuantumCircuit(self.profile.num_orb*2, self.profile.num_orb*2)
        for i in range(self.profile.num_elec//2):
            qc.x(i)
            qc.x(i + self.profile.num_orb)
        self.ansatz.ansatz(qc, self.profile, coeff)
        return qc

    @staticmethod
    def verbose_print(decorator):
        'print output regarding the verbose flag.'
        def wrapper(func):
            def inner_wrapper(self, *args, **kwargs):
                if self.verbose:
                    return decorator(func)(self, *args, **kwargs)
                return func(self, *args, **kwargs)
            return inner_wrapper
        return wrapper

    @staticmethod
    def batch_output(func):
        'Decorator for the batch output.'
        def wrapper(self, *args, **kwargs):
            energy = func(self, *args, **kwargs)
            print(f"Iteration: {self.iteration()}, " +
                  f"Energy: {self.profile.energy_total():12.09f}", end='\r', flush=True)
            return energy
        return wrapper

    @verbose_print(batch_output)
    def batch(self, coeff):
        'Performs the calculation for a given set of coefficients.'
        return self._batch(coeff)

    def _batch(self, coeff):
        qc = self.circuit(coeff)
        energy = self.simulator.measure(qc, self.hamiltonian, self.parallel)
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
            print(f"Iteration: {self.iteration()}, Converged!!         ")
            print(f'Total Energy: {result.energy_total():12.09f}\n')
            print(f'Elapsed time: {elapsed.split(".", maxsplit=1)[0]}')
            return result
        return wrapper

    @verbose_print(general_output)
    def run(self):
        'Performs the VQE calculation.'
        return self._run()

    def _run(self):
        coeff = self.ansatz.generate_coeff(self.profile)
        optimized = self.ansatz.call_optimizer(self.batch, coeff, self.optimizer)
        self.profile.energy_elec = float(optimized.fun)
        self.profile.coeff = optimized.x
        self.profile.circuit = self.circuit(optimized.x)
        # self.profile.spin = self.simulator.measure_spin(self.profile)
        return self.profile

