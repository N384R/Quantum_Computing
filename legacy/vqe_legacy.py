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
from pyscf.gto import Mole
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from jqc.measure.measure import measure
from jqc.measure.hamiltonian import pyscf_luncher, hamiltonian
from ..jqc.vqe.profile import Profile

class VQE:
    'Class for running Variational Quantum Eigensolver (VQE) on a given molecule.'

    def __init__(self, mol: Mole, ansatz = None, simulator = None):
        self.mol: Mole = mol
        self.ansatz = ansatz
        self.simulator = simulator

        self.profile: Profile = Profile(mol)
        self._config = {'iteration': 0, 'verbose': 2, 'parallel': False, 'optimizer': 'powell'}

    @property
    def ansatz(self):
        'The ansatz to be used for the calculation.'
        return self.__ansatz

    @ansatz.setter
    def ansatz(self, ansatz):
        self.__ansatz = ansatz

    @property
    def simulator(self):
        'The simulator to be used for the calculation.'
        return self.__simulator

    @simulator.setter
    def simulator(self, simulator):
        self.__simulator = simulator

    @property
    def haniltonian(self):
        'The Hamiltonian of the molecule.'
        return self.__hamiltonian

    @haniltonian.setter
    def haniltonian(self, haniltonian):
        self.__hamiltonian = haniltonian

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
        self.hamiltonian = hamiltonian(self.profile)
        self._talk('Done')

        self._talk(f'SCF Electronic Energy: {self.profile.energy_elec:18.15f}')
        self._talk(f'SCF Total Energy:      {self.profile.energy_total():18.15f}\n')

        if not self.ansatz:
            raise ValueError('Please provide an ansatz for the calculation.')

    def _initialize_circuit(self, qc):
        for qubit in range(self.profile.num_elec//2):
            qc.x(qubit)
            qc.x(qubit+self.profile.num_orb)

    def _batch(self, coeff):
        self._config['iteration'] += 1
        self._talk(f"\nIteration: {self._config['iteration']}")
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
        energy = measure(qc, self.hamiltonian, self.simulator, self.parallel)
        self._talk('Done')

        self._talk('coeff: ', end='', verb=3)
        self._talk([f'{val:9.6f}' for val in coeff], verb=3)
        self._talk(f'Electronic energy: {energy:18.15f}')
        if self.verbose == 1:
            print(f"Iteration: {self._config['iteration']}, Energy: {energy:18.15f}",
                  end='\r', flush=True)
        return energy

    def run(self):
        'Run calculation'
        start_time = datetime.datetime.now()
        self._init_setup()
        self._talk(f'\nStarting {self.__class__.__name__} Calculation\n')
        self._talk(f'Ansatz: {self.ansatz.__class__.__name__}')
        self._talk(f'Optimizer: {self.optimizer}')

        coeff = self.ansatz.generate_coeff(self.profile)
        optimized = self.ansatz.call_optimizer(self._batch, coeff, self.optimizer)
        self._talk('\n\n!!Successfully Converged!!\n')
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
        if self.verbose == 1:
            print(f"\nOptimized Electronic Energy: {self.profile.energy_elec:18.15f}", flush=True)
            print(f"Optimized Total Energy     : {self.profile.energy_total():18.15f}", flush=True)

        return self.profile

    def _talk(self, line, end='\n', verb=2):
        if verb <= self.verbose:
            print(line, end=end, flush=True)

    def _measure_spin(self, qc):
        shots = 10000
        qc_2 = qc.copy()
        for k in range(self.profile.num_orb * 2):
            qc_2.measure(k, k)

        backend = AerProvider().get_backend('qasm_simulator')
        result = cast(dict[str, float],
                          backend.run(qc_2, shots=shots).result().get_counts())
        self._talk(result)
        spin = 0
        for key, value in result.items():
            for i, orb in enumerate(reversed([*key])):
                if i < 2:
                    spin += value * int(orb)
                else:
                    spin -= value * int(orb)

        return abs(spin/shots/2)
