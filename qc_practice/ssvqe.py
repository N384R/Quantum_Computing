'''
This module contains the implementation of the SSVQE class.
The class is used to calculate the excited state energies of a molecule using the SSVQE algorithm.

The class has the following methods:

__init__: Initializes the SSVQE class with the molecule and the ansatz to be used.
_init_setup: Initializes the active space and the weights to be used in the calculation.
_electron_excitation: Generates the electron excitation operators for the active space.
_swap_init: Swaps the qubits in the circuit to prepare the initial state for the calculation.
_initialize_circuit: Initializes the quantum circuit for the calculation.
_ssvqe_batch: Performs the SSVQE calculation for a given set of coefficients.
run: Runs the SSVQE calculation.
_profiles_update: Updates the profiles with the calculated energies.
'''

from itertools import combinations
import scipy.optimize as opt
from qc_practice import VQE
from qc_practice.profile import Profiles

class SSVQE(VQE):
    '''
    Class for running Subspace-Search Variational Quantum Eigensolver (SSVQE) on a given molecule.
    '''
    def __init__(self, mol, ansatz=None, **kwargs):
        super().__init__(mol, ansatz=ansatz)

        self.active_space = kwargs.get('active_space', None)
        self.weights = kwargs.get('weights', None)
        self.koopmans = kwargs.get('koopmans', False)
        self.verbose = kwargs.get('verbose', 1)

        self.__p = Profiles()
        self.__nspace = None
        self.__transition = None
        self.__trial = 0

    def _init_setup(self):
        super()._init_setup()

        if self.active_space is None:
            self.active_space = [1, 1]

        k = self.active_space
        if ((k[1] > self.profile.num_orb/2) or (k[0] > self.profile.num_elec/2)):
            raise ValueError('Invalid active space. Please check the basis.')

        self.__nspace = k[0] * k[1] * 4
        if not self.koopmans:
            self.__nspace += k[0] * (2 * k[0] - 1) * k[1] * (2 * k[1] - 1)

        if not self.weights:
            self.weights = [1]
            for i in range(1, self.__nspace+1):
                self.weights.append(self.weights[i-1] / 4)

    def _electron_excitation(self, active_space):
        n = self.profile.num_elec//2 - active_space[0]
        occ, vir = active_space
        orbital_indices = list(range(self.profile.num_orb*2))
        alpha_indices = orbital_indices[:self.profile.num_orb]
        beta_indices = orbital_indices[self.profile.num_orb:]
        occ_indices = alpha_indices[:occ] + beta_indices[:occ]
        vir_indices = alpha_indices[occ:vir+1] + beta_indices[occ:vir+1]

        for i in occ_indices:
            for j in vir_indices:
                yield (i+n, j+n)

        if not self.koopmans:
            for i in combinations(occ_indices, 2):
                for j in combinations(vir_indices, 2):
                    k = (i[0] + n, i[1] + n)
                    l = (j[0] + n, j[1] + n)
                    yield (k, l)

    def _swap_init(self, qc, state):
        state = next(self.__transition)
        qc.swap(state[0], state[1])

    def _initialize_circuit(self, qc):
        super()._initialize_circuit(qc)
        if self.profile.state != 0:
            self._swap_init(qc, self.profile.state)

    def _ssvqe_batch(self, coeff):
        self.__trial += 1
        self._talk(f"\nIteration: {self.__trial}")
        self.__transition = self._electron_excitation(self.active_space)
        for i in range(self.__nspace+1):
            self.profile.state = i
            self.__p[i].state = i
            self.verbose, verbose = 0, self.verbose
            self.__p[i].energy_elec = self._batch(coeff)
            self.__p[i].circuit = self.profile.circuit
            self.verbose = verbose
            self._talk(f'State_{i} Energy: {self.__p[i].energy_total():18.15f}')
        weight_en = sum(self.weights[i] * self.__p[i].energy_elec for i in range(self.__nspace+1))

        return weight_en

    def run(self, shots=10000):
        self._init_setup()
        self._talk('\nStarting SSVQE Calculation')
        self.__p.add(self.profile, self.__nspace+1)
        self._shots = shots
        coeff = self.ansatz.generate_coeff(self.profile)
        optimized = opt.minimize(self._ssvqe_batch, coeff, method='SLSQP', tol=1e-6)
        self._talk('\n!!Successfully Converged!!')
        self.profile.coeff = optimized.x
        self.profile = self._profiles_update()
        self._talk('\nFinal Excited State Energies:')
        for i in range(self.__nspace+1):
            self._talk(f'State_{i} Energy: {self.__p[i].energy_total():18.15f}')
        return self.profile

    def _profiles_update(self):
        for i in range(self.__nspace+1):
            self.__p[i].coeff = self.profile.coeff
            self.__p[i].spin = self._measure_spin(self.__p[i].circuit)
        return self.__p
