from pyscf import scf, ao2mo
import scipy.optimize as opt
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from mapper import JordanWignerMapper
from . import Profile

class VQE:
    def __init__(self, mol, verbose=1):
        self.mol = mol
        self.profile = Profile()

        self.verbose = verbose
        self.just_hf = False
        self._shots = None
        self._iterations = 0

        self._talk('Computing Hamiltonian...... ', end='')
        self._initialize_profile()
        self._talk('Done')

        self._talk(f'SCF Electronic Energy: {self.profile.energy_elec:18.15f}')
        self._talk(f'SCF Total Energy:      {self.profile.energy_total():18.15f}\n')

    def _initialize_profile(self):
        rhf = scf.RHF(self.mol)
        rhf.verbose = 0
        rhf.kernel()
        c = rhf.mo_coeff
        hcore = rhf.get_hcore()
        hcore_mo = c.T @ hcore @ c
        self.profile.num_orb = hcore.shape[0]
        self.profile.num_elec = self.mol.nelectron
        N = self.profile.num_orb

        two_elec = self.mol.intor('int2e')
        two_elec_mo = ao2mo.kernel(self.mol, c, two_elec, compact=False)
        two_elec_mo = two_elec_mo.reshape((N, N, N, N))

        self.hamiltonian_pauli = self.hamiltonian(hcore_mo, two_elec_mo)
        self.profile.energy_elec = rhf.energy_elec()[0]
        self.profile.energy_nucl = self.mol.energy_nuc()

    def uccsd_ansatz(self, coeff):
        uccsd_fermion = ''
        N = self.profile.num_orb
        value = iter(coeff)
        for i in range(N//2):
            for j in range(N//2, N):
                val = next(value)
                sign = '+' if val > 0 else '-'
                uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j} ' + '\n'
                sign = '-' if val > 0 else '+'
                uccsd_fermion += f'{sign} {abs(val):f} {j}^ {i} ' + '\n'
                val = next(value)
                sign = '+' if val > 0 else '-'
                uccsd_fermion += f'{sign} {abs(val):f} {i+N}^ {j+N} ' + '\n'
                sign = '-' if val > 0 else '+'
                uccsd_fermion += f'{sign} {abs(val):f} {j+N}^ {i+N} ' + '\n'

        for i in range(N//2):
            for j in range(N//2):
                for k in range(N//2, N):
                    for l in range(N//2, N):
                        val = next(value)
                        sign = '+' if val > 0 else '-'
                        uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j+N}^ {k} {l+N} ' + '\n'
                        sign = '-' if val > 0 else '+'
                        uccsd_fermion += f'{sign} {abs(val):f} {k}^ {l+N}^ {i} {j+N} ' + '\n'

        return JordanWignerMapper(uccsd_fermion)

    def hamiltonian(self, hcore_mo, two_elec_mo):
        second_q = ''
        N = self.profile.num_orb
        for i in range(N):
            coeff = hcore_mo[i, i]
            sign = '+' if coeff > 0 else ''
            second_q += f'{sign} {coeff:.16f} {i}^ {i}' + '\n'
            second_q += f'{sign} {coeff:.16f} {i+N}^ {i+N}' + '\n'

        for i in range(N):
            for j in range(N):
                for k in range(N):
                    for l in range(N):
                        coeff = two_elec_mo[i, j, k, l]/2
                        if coeff < 1e-10:
                            continue
                        sign = '+' if coeff > 0 else ''
                        second_q += f'{sign} {coeff:.16f} {i}^ {k}^ {l} {j}' + '\n'
                        second_q += f'{sign} {coeff:.16f} {i}^ {k+N}^ {l+N} {j}' + '\n'
                        second_q += f'{sign} {coeff:.16f} {i+N}^ {k}^ {l} {j+N}' + '\n'
                        second_q += f'{sign} {coeff:.16f} {i+N}^ {k+N}^ {l+N} {j+N}' + '\n'

        return JordanWignerMapper(second_q)

    def _initialize(self, qc):
        for qubit in range(self.profile.num_elec//2):
            qc.x(qubit)
            qc.x(qubit+self.profile.num_orb)

    def _circuit(self, qc, uccsd_ansatz):
        for p_string, values in uccsd_ansatz.items():
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
        for p_string, values in self.hamiltonian_pauli.items():
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
        self._iterations += 1
        self._talk(f'Iteration: {self._iterations}')
        self._talk('Computing uccsd ansatz..... ', end='')
        uccsd_ansatz = self.uccsd_ansatz(coeff)
        self._talk('Done')

        self._talk('Building quantum circuit... ', end='')
        qc = QuantumCircuit(2*self.profile.num_orb, 2*self.profile.num_orb)
        self._initialize(qc)
        if not self.just_hf:
            self._circuit(qc, uccsd_ansatz)
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
        self._talk('Running VQE\n')
        n = self.profile.num_orb
        self._shots = shots
        coeff = [1e-5] * ((2 * (n//2) **2) + 2 * (n//2 * (n//2 - 1) // 2)**2 + (n//2)**4)
        optimized = opt.minimize(self._batch, coeff, method='Powell')
        self.profile.energy_elec = optimized.fun
        self.profile.coeff = optimized.x
        self._talk(f'Nuclear Repulsion Energy   : {self.profile.energy_nucl:18.15f}')
        self._talk(f'Optimized Electronic Energy: {self.profile.energy_elec:18.15f}\n')
        self._talk(f'Total Energy: {self.profile.energy_total():18.15f}')
        return self.profile.energy_total(), self.profile.coeff

    def _talk(self, line, end='\n', verb=1):
        if verb == self.verbose:
            print(line, end=end)

    def spin(self, qc):
        for k in range(self.profile.num_orb * 2):
            qc.measure(k, k)

        backend = AerProvider().get_backend('qasm_simulator')
        result = backend.run(qc, shots=self._shots).result().get_counts()

        spin = 0
        for key, value in result.items():
            for i, orb in enumerate(reversed([*key])):
                if i < 2:
                    spin += value * int(orb)
                else:
                    spin -= value * int(orb)

        return abs(spin/self._shots/2)
