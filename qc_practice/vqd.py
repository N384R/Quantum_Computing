from qc_practice import VQE
from qc_practice.simulator import Simulator
from qc_practice.simulator import QASM
from .profile import Profiles

class VQD(VQE):
    '''
    Variational Quantum Deflation (VQD).

    Example:
    >>> from pyscf import gto
    >>> from qc_practice.ansatz import UCCSD
    >>> from qc_practice.simulator import QASM

    '''
    def __init__(self, mol, ansatz = None, simulator: Simulator = QASM()):
        super().__init__(mol, ansatz = ansatz, simulator = simulator)
        self._profiles = None  #type: ignore
        self.nstates = 2

    @property
    def nstates(self):
        'The number of states for the calculation.'
        return self._config['nstates']

    @nstates.setter
    def nstates(self, nstates):
        self._config['nstates'] = nstates

    def batch(self, coeff):
        'Performs a batch calculation.'
        energy = super().batch(coeff)
        beta = 3.0
        for i in range(self.profile.state):
            state1 = self.profile
            state2 = self._profiles[i]
            overlap_sq = self.simulator.swap_test(state1, state2, self.ansatz)
            energy += beta * overlap_sq
        return energy

    @staticmethod
    def general_output(func):
        'Decorator for the normal output.'
        def wrapper(self, *args, **kwargs):
            print(f'\nStarting {self.__class__.__name__} Calculation\n')
            print(f'Ansatz: {self.ansatz.__class__.__name__}')
            print(f'Simulator: {self.simulator.__class__.__name__}')
            print(f'Optimizer: {self.optimizer}')
            print(f'nstates: {self.nstates}\n')
            result = func(self, *args, **kwargs)
            print('Final State Energies:')
            print(self.profile.energy_total())
            return result
        return wrapper

    @general_output
    def run(self):
        'Performs the VQD calculation.'
        self._profiles: Profiles = Profiles(self.profile, self.nstates)
        for i in range(self.nstates):
            self.profile.state = i
            self._config['iteration'] = 0
            super()._run()
            self._profiles.update(self.profile)
        self.profile = self._profiles  #type: ignore
        return self.profile
