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

import datetime
from typing import cast
from itertools import combinations
import scipy.optimize as opt
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import VQE
from .profile import Profiles

class SSVQE(VQE):
    '''
    Class for running Subspace-Search Variational Quantum Eigensolver (SSVQE) on a given molecule.
    '''
    def __init__(self, mol, ansatz = None):
        super().__init__(mol, ansatz = ansatz)

        self._p = Profiles()
        self._config['trial'] = 0
        self._config['active_space'] = [1, 1]
        self._config['koopmans'] = False
        self._config['weights'] = None
        self._transition = iter(())

    @property
    def active_space(self):
        'The active space for the calculation.'
        return self._config['active_space']

    @active_space.setter
    def active_space(self, active_space):
        self._config['active_space'] = active_space

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

    def _init_setup(self):
        super()._init_setup()

        k = self.active_space
        no = self.profile.num_orb
        ne = self.profile.num_elec
        if ((k[1] > no - ne//2) or (k[0] > ne//2)):
            raise ValueError('Invalid active space. Please check the basis.')

        self._config['nspace'] = k[0] * k[1] * 4
        if not self.koopmans:
            self._config['nspace'] += k[0] * (2 * k[0] - 1) * k[1] * (2 * k[1] - 1)

        if not self.weights:
            self.weights = [self._config['nspace']+5]
            for i in range(1, self._config['nspace'] + 1):
                self.weights.append(self._config['nspace'] - i + 1)

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

    def _swap_init(self, qc):
        state = next(self._transition)
        print(state)
        qc.swap(state[0], state[1])

    def _initialize_circuit(self, qc):
        super()._initialize_circuit(qc)
        if self._config['state'] != 0:
            self._swap_init(qc)

    def _excite_batch(self, coeff):
        self._config['trial'] += 1
        self._talk(f"\nIteration: {self._config['trial']}")
        self._transition = self._electron_excitation(self.active_space)
        for i in range(self._config['nspace']+1):
            self._config['state'] = i
            self.verbose, verbose = 0, self.verbose
            self._p[i].energy_elec = self._batch(coeff)
            self._p[i].transition = self._transition
            self._p[i].circuit = self.profile.circuit
            self._p[i].state = i
            self.verbose = verbose
            self._talk(f'State_{i} Energy: {self._p[i].energy_total():18.15f}')
        weight_en = sum(self.weights[i] * self._p[i].energy_elec
                        for i in range(self._config['nspace']+1))

        return weight_en

    def run(self, shots=10000):
        start_time = datetime.datetime.now()
        self._init_setup()
        self._talk(f'\nStarting {self.__class__.__name__} Calculation\n')
        self._talk(f'Ansatz: {self.ansatz.__class__.__name__}')
        self._talk(f'Optimizer: {self.optimizer}')
        self._talk(f'Active Space: {self.active_space}')
        self._talk(f'Koopmans Condition: {self.koopmans}')
        self._talk(f'Weights: {self.weights}')

        self._p.add(self.profile, self._config['nspace']+1)
        self._config['shots'] = shots
        coeff = self.ansatz.generate_coeff(self.profile)
        optimized = opt.minimize(self._excite_batch, coeff, method=self.optimizer)
        self._talk('\n!!Successfully Converged!!')
        self.profile.coeff = optimized.x
        self.profile = self._profiles_update() # type: ignore
        self._talk('\nFinal Excited State Energies:')
        for i in range(self._config['nspace']+1):
            self._talk(f'State_{i} Energy: {self._p[i].energy_total():18.15f}')

        elapsed_time = str(datetime.datetime.now() - start_time)
        self._talk(f'\nElapsed Time    : {elapsed_time.split(".", maxsplit=1)[0]}')

        return self.profile

    def _profiles_update(self):
        for i in range(self._config['nspace']+1):
            self._p[i].coeff = self.profile.coeff
            self._p[i].spin = self._measure_spin(self._p[i].circuit)
        return self._p

    def _calculate_overlap(self):
        no = self.profile.num_orb
        qc = [QuantumCircuit(no*2, no*2) for _ in range(2)]
        ancila = QuantumCircuit(1, 1)
        for i in range(2):
            self._config['state'] = i
            self._initialize_circuit(qc[i])
            self.ansatz.ansatz(qc[i], self._p[i], self._p[i].coeff)
        qc_total = ancila.tensor(qc[0]).tensor(qc[1])

        qc_total.h(0)
        for i in range(1, no*2+1):
            qc_total.ccx(i+no*2, i, 0)
        qc_total.h(0)

        backend = AerProvider().get_backend('qasm_simulator')
        result = cast(dict[str, float],
                        backend.run(qc_total, shots=self._config['shots']).result().get_counts())
        print(result)
