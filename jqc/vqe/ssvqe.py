import datetime
from itertools import combinations
from qiskit import QuantumCircuit
from jqc.vqe.vqe import VQE
from jqc.ansatz import Ansatz
from jqc.ansatz import fUCCSD
from jqc.simulator import Simulator
from jqc.simulator import StateVector
from jqc.vqe.profile import Profiles


class SSVQE(VQE):
    '''
    Subspace-Search Variational Quantum Eigensolver (SSVQE).
    '''
    verbose_print = VQE.verbose_print

    def __init__(self, mol, ansatz: Ansatz = fUCCSD(), simulator: Simulator = StateVector()):
        super().__init__(mol, ansatz=ansatz, simulator=simulator)
        self.profiles = None  # type: ignore
        self.koopmans = False
        self.active_space = [1, 1]

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
                raise ValueError(
                    'Invalid active space. Please check the basis.')
            self._config['nstates'] = len(self._configuration())
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
        if 'weights' not in self._config:
            self._config['weights'] = [self.nstates * 10] + \
                list(range(1, self.nstates))[::-1]
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
        nf = ne//2 + vir
        orb_indices = list(range(no * 2))
        alpha = orb_indices[ni:nf]
        beta = orb_indices[no + ni:no + nf]
        idx = alpha + beta
        configuration = [idx for i in range(ni) for idx in (i, i + no)]
        return [configuration + list(k) for k in combinations(idx, occ * 2)]

    def circuit(self, coeff):
        'Builds the quantum circuit for the calculation.'
        qc = QuantumCircuit(self.profile.num_orb*2, self.profile.num_orb*2)
        for idx in self._configuration()[self.profile.state]:
            qc.x(idx)
        self.ansatz.ansatz(qc, self.profile, coeff)
        return qc

    @staticmethod
    def state_output(func):
        'Decorator for the batch output.'

        def wrapper(self, *args, **kwargs):
            energy = func(self, *args, **kwargs)
            print(
                f"State {self.profile.state}: {self.profile.energy_total:12.09f}")
            return energy
        return wrapper

    @verbose_print(state_output)
    def _batch(self, coeff):
        'Performs the calculation for a given set of coefficients.'
        return super()._batch(coeff)

    @staticmethod
    def batch_output(func):
        'Decorator for the batch output.'

        def wrapper(self, *args, **kwargs):
            print(f"Iteration: {self.iteration()}")
            energy = func(self, *args, **kwargs)
            print(f"\x1b[{self.nstates+2}A")
            return energy
        return wrapper

    @verbose_print(batch_output)
    def batch(self, coeff):
        for state in range(self.nstates):
            self.profile.state = state
            self._batch(coeff)
            self.profiles.update(self.profile)
        cost_func = sum(self.weights[i] * self.profiles[i].energy_elec
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
            print(f"\x1b[{self.nstates+1}B")
            print('!!Successfully Converged!!\n')
            print('Final State Energies:')
            for i, state in enumerate(self.profiles):
                print(f'State {i}: {state.energy_total:12.09f}')
            print(f'\nElapsed time: {elapsed.split(".", maxsplit=1)[0]}')
            del self.profile
            return result
        return wrapper

    @verbose_print(general_output)
    def run(self):
        'Performs the SSVQE calculation.'
        self.profiles: Profiles = Profiles(self.profile, self.nstates)
        super()._run()
        return self

    def transition_matrix(self, operator):
        'Calculates the transition matrix.'
        def transition_matrix_element(state1, state2):
            def measure_value(mode):
                qc = superposition(state1, state2, mode)
                self.ansatz.ansatz(qc, state1, state1.coeff)
                result = self.simulator.measure(qc, operator, self.parallel)
                result -= self.simulator.measure(state1.circuit,
                                                 operator, self.parallel) / 2
                result -= self.simulator.measure(state2.circuit,
                                                 operator, self.parallel) / 2
                return result
            real = measure_value('real')
            imag = measure_value('imag')
            return real + 1j * imag

        matrix = [[transition_matrix_element(state1, state2)
                   for state1 in self.profiles]
                  for state2 in self.profiles]
        return matrix

# def superposition(state1, state2, mode):
#     'Creates a superposition of two states.'
#     no = state1.num_orb
#     c1 = state1.configuration
#     c2 = state2.configuration
#     qc = QuantumCircuit(no*2, no*2)
#     if c1[0] != c2[0]:
#         qc.h(c1[0])
#         if mode == 'imag':
#             qc.sdg(c1[0])
#         qc.cx(c1[0], c2[0])
#         qc.x(c2[0])
#     else:
#         qc.x(c1[0])
#     if c1[1] != c2[1]:
#         qc.h(c1[1])
#         if mode == 'imag':
#             qc.sdg(c1[1])
#         qc.cx(c1[1], c2[1])
#         qc.x(c2[1])
#     else:
#         qc.x(c1[1])
#     return qc
