from pyscf import scf, ao2mo
import scipy.optimize as opt
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import JordanWignerMapper

class VQE:
    def __init__(self, mol, verbose=1):
        self.mol = mol
        rhf = scf.RHF(self.mol)
        rhf.verbose = 0
        rhf.kernel()
        c = rhf.mo_coeff
        hcore = rhf.get_hcore()
        hcore_mo = c.T @ hcore @ c
        self._num = hcore.shape[0]
        self._num_elec = self.mol.nelectron
        self.nuclear_repulsion = self.mol.energy_nuc()
        two_elec = self.mol.intor('int2e')
        # TWO_ELEC_MO = np.einsum('pi,qj,rk,sl,ijkl->pqrs', C, C, C, C, mol.intor('int2e'))
        two_elec_mo = ao2mo.kernel(self.mol, c, two_elec, compact=False)
        two_elec_mo = two_elec_mo.reshape((self._num, self._num, self._num, self._num))

        self.verbose = verbose
        self._shots = None
        self._iterations = 0

        self._talk('Computing Hamiltonian...... ', end='')
        self.hamiltonian_pauli = self.hamiltonian(hcore_mo, two_elec_mo)
        self._talk('Done')
        total_energy = rhf.energy_elec()[0] + self.nuclear_repulsion
        self._talk(f'SCF Electronic Energy: {rhf.energy_elec()[0]:18.15f}')
        self._talk(f'SCF Total Energy:      {total_energy:18.15f}\n')

    def uccsd_ansatz(self, coeff):
        uccsd_fermion = ''
        N = self._num
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
        N = self._num
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
        for qubit in range(self._num_elec//2):
            qc.x(qubit)
            qc.x(qubit+self._num)

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
        qc = QuantumCircuit(2*self._num, 2*self._num)
        self._initialize(qc)
        self._circuit(qc, uccsd_ansatz)
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
        n = self._num
        self._shots = shots
        coeff = [1e-5] * ((2 * (n//2) **2) + 2 * (n//2 * (n//2 - 1) // 2)**2 + (n//2)**4)
        optimized_energy = opt.minimize(self._batch, coeff, method='Powell')
        total_energy = optimized_energy.fun + self.nuclear_repulsion
        self._talk(f'Nuclear Repulsion Energy   : {self.nuclear_repulsion:18.15f}')
        self._talk(f'Optimized Electronic Energy: {optimized_energy.fun:18.15f}\n')
        self._talk(f'Total Energy: {total_energy:18.15f}')
        return total_energy, optimized_energy.x

    def run_hf(self, shots=10000):
        self._shots = shots
        qc = QuantumCircuit(2*self._num, 2*self._num)
        for qubit in range(self._num_elec//2):
            qc.x(qubit)
            qc.x(qubit+self._num)

        self._talk('Measuring energy... ', end='')
        energy = self._measure(qc)
        self._talk('Done')

        total_energy = energy + self.nuclear_repulsion
        self._talk(f'Optimized Electronic Energy: {total_energy:18.15f}', 1)
        return total_energy

    def _talk(self, line, end='\n', verb=1):
        if verb == self.verbose:
            print(line, end=end)
