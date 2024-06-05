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

from pyscf import scf, ao2mo
import scipy.optimize as opt
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice.mapper.jordan_wigner import JordanWignerMapper
from .profile import Profile

class VQE:
    'Class for running Variational Quantum Eigensolver (VQE) on a given molecule.'

    def __init__(self, mol, ansatz=None, verbose=1):
        self.mol = mol
        self.ansatz = ansatz
        self.profile = Profile()
        self.verbose = verbose
        self.just_hf = False
        self._shots = None

        self.__iteration = 0
        self.__hamiltonian_pauli = None

    def _init_setup(self):
        self._talk('Computing Hamiltonian...... ', end='')
        self._initialize_profile()
        self._talk('Done')

        self._talk(f'SCF Electronic Energy: {self.profile.energy_elec:18.15f}')
        self._talk(f'SCF Total Energy:      {self.profile.energy_total():18.15f}\n')

        if not self.ansatz:
            raise ValueError('Please provide an ansatz for the calculation.')

    def _initialize_profile(self):
        rhf = scf.RHF(self.mol)
        rhf.verbose = 0
        rhf.kernel()
        c = rhf.mo_coeff
        hcore = rhf.get_hcore()
        hcore_mo = c.T @ hcore @ c
        self.profile.num_orb = hcore.shape[0]
        self.profile.num_elec = self.mol.nelectron
        n = self.profile.num_orb

        two_elec = self.mol.intor('int2e')
        two_elec_mo = ao2mo.kernel(self.mol, c, two_elec, compact=False)
        two_elec_mo = two_elec_mo.reshape((n, n, n, n))

        self.__hamiltonian_pauli = self._hamiltonian(hcore_mo, two_elec_mo)
        self.profile.energy_elec = rhf.energy_elec()[0]
        self.profile.energy_nucl = self.mol.energy_nuc()

    def _hamiltonian(self, hcore_mo, two_elec_mo):
        second_q = ''
        n = self.profile.num_orb
        for i in range(n):
            coeff = hcore_mo[i, i]
            sign = '+' if coeff > 0 else ''
            second_q += f'{sign} {coeff:.16f} {i}^ {i}' + '\n'
            second_q += f'{sign} {coeff:.16f} {i+n}^ {i+n}' + '\n'

        for i in range(n):
            for j in range(n):
                for k in range(n):
                    for l in range(n):
                        coeff = two_elec_mo[i, j, k, l]/2
                        if coeff < 1e-10:
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

    def _circuit(self, qc, ansatz):
        for p_string, values in ansatz.items():
            chk = []
            for idx, p in p_string.items():
                if p.symbol == 'X':
                    qc.h(idx)

                elif p.symbol == 'Y':
                    qc.s(idx)
                    qc.h(idx)

                if p.symbol != 'I':
                    chk += [idx]

            for i in chk:
                if i != max(chk):
                    qc.cx(i, i+1)

            qc.rz(values.imag, max(chk))

            for i in reversed(chk):
                if i != max(chk):
                    qc.cx(i, i+1)

            for idx, p in p_string.items():
                if p.symbol == 'X':
                    qc.h(idx)

                elif p.symbol == 'Y':
                    qc.h(idx)
                    qc.sdg(idx)

    def _measure(self, qc):
        energy = 0.
        for p_string, values in self.__hamiltonian_pauli.items():
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
                energy += values.real
                self._talk(f'Expectation: {values.real:18.15f} {p_string}', verb=2)
                continue

            backend = AerProvider().get_backend('qasm_simulator')
            result = backend.run(qc_2, shots=self._shots).result().get_counts()

            counts = 0
            for key, value in result.items():
                counts += (-1)**sum(int(k) for k in key) * value

            expectation = counts/self._shots * values.real

            self._talk(f'Expectation: {expectation:18.15f} {p_string}', verb=2)
            energy += expectation
        return energy

    def _batch(self, coeff):
        self.__iteration += 1
        self._talk(f"Iteration: {self.__iteration}")
        self._talk('Computing uccsd ansatz..... ', end='')
        ansatz = self.ansatz.ansatz(self.profile, coeff)
        self._talk('Done')

        self._talk('Building quantum circuit... ', end='')
        qc = QuantumCircuit(2*self.profile.num_orb, 2*self.profile.num_orb)
        self._initialize_circuit(qc)
        if not self.just_hf:
            self._circuit(qc, ansatz)
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
        'Run the VQE optimization'
        self._init_setup()
        self._talk('Starting VQE Optimization... ')

        self._shots = shots
        coeff = self.ansatz.generate_coeff(self.profile)
        optimized = opt.minimize(self._batch, coeff, method='COBYLA', tol=1e-6)
        self._talk('\n!!Successfully Converged!!\n')
        self.profile.energy_elec = optimized.fun
        self.profile.coeff = optimized.x

        self.profile.spin = self._measure_spin(self.profile.circuit)

        self._talk(f'Nuclear Repulsion Energy   : {self.profile.energy_nucl:18.15f}')
        self._talk(f'Optimized Electronic Energy: {self.profile.energy_elec:18.15f}\n')
        self._talk(f'Total Energy        : {self.profile.energy_total():18.15f}')
        self._talk(f'Coefficients        : {self.profile.coeff}')
        self._talk(f'Multiplicity (2S+1) : {self.profile.spin:.02f}')

        return self.profile

    def _talk(self, line, end='\n', verb=1):
        if verb == self.verbose:
            print(line, end=end)

    def _measure_spin(self, qc):
        for k in range(self.profile.num_orb * 2):
            qc.measure(k, k)

        backend = AerProvider().get_backend('qasm_simulator')
        result = backend.run(qc, shots=self._shots).result().get_counts()
        # print(result)
        spin = 0
        for key, value in result.items():
            for i, orb in enumerate(reversed([*key])):
                if i < 2:
                    spin += value * int(orb)
                else:
                    spin -= value * int(orb)

        return abs(spin/self._shots/2)
