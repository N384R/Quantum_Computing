'''
This module contains the implementation of the SSVQE class. 
The class is initialized with the molecule and the ansatz to be used.

The class has the following methods:

__init__: Initializes the VQE class with the molecule and the ansatz to be used.
_initialize_profile: Initializes the profile with the Hamiltonian and the electronic energy.
_hamiltonian: Generates the Hamiltonian in the Pauli basis.
_initialize_circuit: Initializes the quantum circuit for the calculation.
_circuit: Builds the quantum circuit for the calculation.
_measure: Measures the energy of the quantum circuit.
_batch: Performs the VQE calculation for a given set of coefficients.
run: Runs the VQE optimization.
_talk: Prints the output based on the verbosity level.
_measure_spin: Measures the spin of the quantum circuit.
'''

import datetime
from typing import cast
from itertools import product
from multiprocessing import Pool
import numpy as np
from pyscf import ao2mo, scf
from pyscf.gto import Mole
import scipy.optimize as opt
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice.mapper.jordan_wigner import JordanWignerMapper
from .profile import Profile

class VQE:
    'Class for running Variational Quantum Eigensolver (VQE) on a given molecule.'

    def __init__(self, mol: Mole, ansatz = None):
        self.mol: Mole = mol
        self.ansatz = ansatz
        self._hamiltonian_pauli = {}
        self.profile: Profile = Profile()
        self._config = {'iteration': 0, 'verbose': 1, 'parallel': False, 'optimizer': 'powell'}

    @property
    def ansatz(self):
        'The ansatz to be used for the calculation.'
        return self.__ansatz

    @ansatz.setter
    def ansatz(self, ansatz):
        self.__ansatz = ansatz

    @property
    def verbose(self):
        'The verbosity level of the calculation.'
        return self._config['verbose']

    @verbose.setter
    def verbose(self, verbose):
        self._config['verbose'] = verbose

    @property
    def parallel(self):
        'The parallelization flag for the calculation.'
        return self._config['parallel']

    @parallel.setter
    def parallel(self, parallel):
        if not isinstance(parallel, bool):
            raise ValueError("Parallel flag must be a boolean.")
        self._config['parallel'] = parallel

    @property
    def optimizer(self):
        'The optimizer to be used for the calculation.'
        return self._config['optimizer']

    @optimizer.setter
    def optimizer(self, optimizer):
        self._config['optimizer'] = optimizer

    def _init_setup(self):
        self._talk('Computing Hamiltonian...... ', end='')
        self._initialize_profile()
        self._talk('Done')

        self._talk(f'SCF Electronic Energy: {self.profile.energy_elec:18.15f}')
        self._talk(f'SCF Total Energy:      {self.profile.energy_total():18.15f}\n')

        if not self.ansatz:
            raise ValueError('Please provide an ansatz for the calculation.')

    def _initialize_profile(self):
        mf = scf.RHF(self.mol)
        mf.verbose = 0
        mf.kernel()

        c = np.asarray(mf.mo_coeff)
        hcore = mf.get_hcore()
        hcore_mo = c.T @ hcore @ c
        self.profile.num_elec = self.mol.nelectron
        self.profile.num_orb = hcore.shape[0]
        n = self.profile.num_orb

        two_elec = self.mol.intor('int2e')
        two_elec_mo = np.asarray(
            ao2mo.kernel(self.mol, c, two_elec, compact=False)
            ).reshape((n, n, n, n))

        self._hamiltonian_pauli = self._hamiltonian(hcore_mo, two_elec_mo)
        self.profile.energy_elec = mf.energy_elec()[0]
        self.profile.energy_nucl = self.mol.energy_nuc()

    def _hamiltonian(self, hcore_mo, two_elec_mo):
        second_q = ''
        n = self.profile.num_orb
        for i, j in product(range(n), repeat=2):
            if abs(coeff := hcore_mo[i, j]) < 1e-10:
                continue
            sign = '+' if coeff > 0 else ''
            second_q += f'{sign} {coeff:.16f} {i}^ {j}' + '\n'
            second_q += f'{sign} {coeff:.16f} {i+n}^ {j+n}' + '\n'

        for i, j, k, l in product(range(n), repeat=4):
            if abs(coeff := two_elec_mo[i, j, k, l]/2) < 1e-10:
                continue
            sign = '+' if coeff > 0 else ''
            second_q += f'{sign} {coeff:.16f} {i}^ {k}^ {l} {j}' + '\n'
            second_q += f'{sign} {coeff:.16f} {i}^ {k+n}^ {l+n} {j}' + '\n'
            second_q += f'{sign} {coeff:.16f} {i+n}^ {k}^ {l} {j+n}' + '\n'
            second_q += f'{sign} {coeff:.16f} {i+n}^ {k+n}^ {l+n} {j+n}' + '\n'

        return JordanWignerMapper(second_q)

    def _initialize_circuit(self, qc):
        for qubit in range(self.profile.num_elec//2):
            qc.x(qubit)
            qc.x(qubit+self.profile.num_orb)

    @staticmethod
    def _single_measure(args: tuple[QuantumCircuit, dict, complex, int]):
        qc, p_string, values, shots = args
        qc_2 = qc.copy()
        for idx, p in p_string.items():
            if p.symbol == 'X':
                qc_2.h(idx)

            elif p.symbol == 'Y':
                qc_2.sdg(idx)
                qc_2.h(idx)

        for idx, p in p_string.items():
            if p.symbol == 'I':
                continue

            qc_2.measure(idx, idx)

        if all(p.symbol == 'I' for p in p_string.values()):
            return values.real

        backend = AerProvider().get_backend('qasm_simulator')
        result = cast(dict[str, float],
                        backend.run(qc_2, shots=shots).result().get_counts())

        counts = 0
        for key, value in result.items():
            counts += (-1)**sum(int(k) for k in key) * value

        expectation = counts / shots * values.real
        return expectation

    def _measure(self, qc):
        tasks = [(qc, p_string, values, self._config['shots'])
                 for p_string, values in self._hamiltonian_pauli.items()]

        if self.parallel:
            with Pool(processes=4) as pool:
                results = pool.map(self._single_measure, tasks)
        else:
            results = [self._single_measure(task) for task in tasks]

        energy = sum(results)
        return energy

    def _batch(self, coeff):
        self._config['iteration'] += 1
        self._talk(f"Iteration: {self._config['iteration']}")
        self._talk('Computing uccsd ansatz..... ', end='')
        self._talk('Done')

        self._talk('Building quantum circuit... ', end='')
        qc = QuantumCircuit(2*self.profile.num_orb, 2*self.profile.num_orb)
        self._initialize_circuit(qc)
        if self.ansatz:
            self.ansatz.ansatz(qc, self.profile, coeff)
        self.profile.circuit = qc
        self._talk('Done')

        self._talk('Measuring energy........... ', end='')
        energy = self._measure(qc)
        self._talk('Done')

        self._talk('coeff: ', end='', verb=2)
        self._talk([f'{val:9.6f}' for val in coeff], verb=2)
        self._talk(f'Electronic energy: {energy:18.15f}\n')
        return energy

    def run(self, shots=10000):
        'Run calculation'
        start_time = datetime.datetime.now()
        self._init_setup()
        self._talk(f'\nStarting {self.__class__.__name__} Calculation\n')
        self._talk(f'Ansatz: {self.ansatz.__class__.__name__}')
        self._talk(f'Optimizer: {self.optimizer}')

        self._config['shots'] = shots
        coeff = self.ansatz.generate_coeff(self.profile)
        optimized = opt.minimize(self._batch, coeff, method=self.optimizer)
        self._talk('\n!!Successfully Converged!!\n')
        self.profile.energy_elec = optimized.fun
        self.profile.coeff = optimized.x

        self.profile.spin = self._measure_spin(self.profile.circuit)
        elapsed_time = str(datetime.datetime.now() - start_time)

        self._talk(f'Nuclear Repulsion Energy   : {self.profile.energy_nucl:18.15f}')
        self._talk(f'Optimized Electronic Energy: {self.profile.energy_elec:18.15f}\n')
        self._talk(f'Total Energy        : {self.profile.energy_total():18.15f}')
        self._talk(f'Coefficients        : {self.profile.coeff}')
        self._talk(f'Multiplicity (2S+1) : {self.profile.spin:.02f}\n')
        self._talk(f'Elapsed Time        : {elapsed_time.split(".", maxsplit=1)[0]}')

        return self.profile

    def _talk(self, line, end='\n', verb=1):
        if verb == self.verbose:
            print(line, end=end, flush=True)

    def _measure_spin(self, qc):
        qc_2 = qc.copy()
        for k in range(self.profile.num_orb * 2):
            qc_2.measure(k, k)

        backend = AerProvider().get_backend('qasm_simulator')
        result = cast(dict[str, float],
                          backend.run(qc_2, shots=self._config['shots']).result().get_counts())
        print(result)
        spin = 0
        for key, value in result.items():
            for i, orb in enumerate(reversed([*key])):
                if i < 2:
                    spin += value * int(orb)
                else:
                    spin -= value * int(orb)

        return abs(spin/self._config['shots']/2)
