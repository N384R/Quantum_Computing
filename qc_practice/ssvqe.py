import datetime
from qiskit import QuantumCircuit
from qc_practice import VQE
from qc_practice.ansatz import Ansatz
from qc_practice.ansatz import UpCCGSD
from qc_practice.simulator import Simulator
from qc_practice.simulator import QASM
from qc_practice.profile import Profiles

class SSVQE(VQE):
    '''
    Subspace-Search Variational Quantum Eigensolver (SSVQE).
    '''
    def __init__(self, mol, ansatz: Ansatz = UpCCGSD(), simulator: Simulator=QASM()):
        super().__init__(mol, ansatz = ansatz, simulator = simulator)
        self._profiles = None  #type: ignore
        self.koopmans = False
        self.active_space = [1, 1]
        self.weights = [self.nstates + 5] + list(reversed(range(1, self.nstates)))
        self.configuration_time = 0

    @property
    def active_space(self):
        'The number of spaces for the calculation.'
        return self._config['active_space']

    @active_space.setter
    def active_space(self, active_space):
        self._config['active_space'] = active_space

    @property
    def nstates(self):
        'The number of states for the calculation.'
        k = self.active_space
        if 'nstates' not in self._config:
            no = self.profile.num_orb
            ne = self.profile.num_elec
            if ((k[1] > no - ne//2) or (k[0] > ne//2)):
                raise ValueError('Invalid active space. Please check the basis.')
        self._config['nstates'] = 1 + k[0] * k[1] * 2
        if not self.koopmans:
            self._config['nstates'] += k[0] * (2 * k[0] - 1) * k[1] * (2 * k[1] - 1)
        return self._config['nstates']

    @property
    def koopmans(self):
        'The Koopmans condition for the calculation.'
        return self._config['koopmans']

    @koopmans.setter
    def koopmans(self, koopmans):
        self._config['koopmans'] = koopmans

    @property
    def weights(self):
        'The weights for the calculation.'
        return self._config['weights']

    @weights.setter
    def weights(self, weights):
        self._config['weights'] = weights

    def _configuration(self):
        'Initializes the configuration.'
        occ, vir = self.active_space
        ne = self.profile.num_elec
        no = self.profile.num_orb
        ni = ne//2 - occ
        nf = ni + vir + 1
        orb_indices = list(range(no * 2))
        alpha = orb_indices[ni:nf]
        beta = orb_indices[no + ni:no + nf]
        return [(i, j) for i in alpha for j in beta]

    def circuit(self, coeff):
        'Builds the quantum circuit for the calculation.'
        qc = QuantumCircuit(self.profile.num_orb*2, self.profile.num_orb*2)
        i, j = self._configuration()[self.profile.state]
        qc.x(i)
        qc.x(j)
        self.ansatz.ansatz(qc, self.profile, coeff)
        return qc

    @staticmethod
    def state_output(func):
        'Decorator for the batch output.'
        def wrapper(self, *args, **kwargs):
            energy = func(self, *args, **kwargs)
            print(f"State {self.profile.state}: {self.profile.energy_total():12.09f}")
            return energy
        return wrapper

    @state_output
    def state_batch(self, coeff):
        'Performs the calculation for a given set of coefficients.'
        return self._batch(coeff)

    @staticmethod
    def batch_output(func):
        'Decorator for the batch output.'
        def wrapper(self, *args, **kwargs):
            print(f"Iteration: {self.iteration()}")
            energy = func(self, *args, **kwargs)
            print(f"\x1b[{self.nstates+2}A")
            return energy
        return wrapper

    @batch_output
    def batch(self, coeff):
        for state in range(self.nstates):
            self.profile.state = state
            self.state_batch(coeff)
            self._profiles.update(self.profile)
        cost_func = sum(self.weights[i] * self._profiles[i].energy_elec
                        for i in range(self.nstates))
        return cost_func

    @staticmethod
    def general_output(func):
        'Decorator for the normal output.'
        def wrapper(self, *args, **kwargs):
            print(f'\nStarting {self.__class__.__name__} Calculation\n')
            print(f'Ansatz: {self.ansatz.__class__.__name__}')
            print(f'Simulator: {self.simulator.__class__.__name__}')
            print(f'Optimizer: {self.optimizer}')
            print(f'Koopmans: {self.koopmans}')
            print(f'nstates: {self.nstates}')
            print(f'Active Space: {self.active_space}')
            print(f'weights: {self.weights}\n')
            start = datetime.datetime.now()
            result = func(self, *args, **kwargs)
            elapsed = str(datetime.datetime.now() - start)
            print('\n'*(self.nstates+1))
            print('!!Successfully Converged!!\n')
            print('Final State Energies:')
            for i, state in enumerate(self.profile):
                print(f'State {i}: {state.energy_total():12.09f}')
            print(f'\nElapsed time: {elapsed.split(".", maxsplit=1)[0]}')
            return result
        return wrapper

    @general_output
    def run(self):
        'Performs the SSVQE calculation.'
        self._profiles: Profiles = Profiles(self.profile, self.nstates)
        super()._run()
        self.profile = self._profiles  #type: ignore
        return self.profile
