import datetime
from qc_practice import VQE
from qc_practice.ansatz import Ansatz
from qc_practice.ansatz import UpCCGSD
from qc_practice.simulator import Simulator
from qc_practice.simulator import QASM
from qc_practice.profile import Profiles

class VQD(VQE):
    '''
    Variational Quantum Deflation (VQD).

    Example:
    >>> from pyscf import gto
    >>> from qc_practice import VQD
    >>> from qc_practice.ansatz import UCCSD
    >>> from qc_practice.simulator import QAS
    >>> mol = gto.M(atom = 'H 0 0 0; H 0 0 0.7', basis = 'sto-3g')
    >>> vqd = VQD(mol, ansatz = UCCSD(), simulator = QASM())
    >>> result = vqd.run()
    '''
    def __init__(self, mol, ansatz: Ansatz = UpCCGSD(), simulator: Simulator = QASM()):
        super().__init__(mol, ansatz = ansatz, simulator = simulator)
        self._profiles = None  #type: ignore
        self.nstates = 2
        self.beta = 3.0

    @property
    def nstates(self):
        'The number of states for the calculation.'
        return self._config['nstates']

    @nstates.setter
    def nstates(self, nstates: int):
        self._config['nstates'] = nstates

    @property
    def beta(self):
        'The beta parameter for the calculation.'
        return self._config['beta']

    @beta.setter
    def beta(self, beta: float):
        self._config['beta'] = beta

    def batch(self, coeff):
        'Performs a batch calculation.'
        energy = super().batch(coeff)
        for i in range(self.profile.state):
            state1 = self.profile
            state2 = self._profiles[i]
            overlap_sq = self.simulator.get_overlap(state1, state2)
            energy += self.beta * overlap_sq
        return energy

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
        return super()._run()

    @staticmethod
    def general_output(func):
        'Decorator for the normal output.'
        def wrapper(self, *args, **kwargs):
            print(f'\nStarting {self.__class__.__name__} Calculation\n')
            print(f'Ansatz: {self.ansatz.__class__.__name__}')
            print(f'Simulator: {self.simulator.__class__.__name__}')
            print(f'Optimizer: {self.optimizer}')
            print(f'nstates: {self.nstates}\n')
            start = datetime.datetime.now()
            result = func(self, *args, **kwargs)
            elapsed = str(datetime.datetime.now() - start)
            print('Final State Energies:')
            for i, state in enumerate(self.profile):
                print(f'State {i}: {state.energy_total():12.09f}')
            print(f'\nElapsed time: {elapsed.split(".", maxsplit=1)[0]}')
            return result
        return wrapper

    @general_output
    def run(self):
        'Performs the VQD calculation.'
        self._profiles: Profiles = Profiles(self.profile, self.nstates)
        for i in range(self.nstates):
            self.profile.state = i
            self._config['iteration'] = 0
            self._run()
            self._profiles.update(self.profile)
        self.profile = self._profiles  #type: ignore
        return self.profile
