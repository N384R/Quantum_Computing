from pyscf import scf, ao2mo
import scipy.optimize as opt
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import JordanWignerMapper

class VQE:
    def __init__(self, mol):
        self.mol = mol
        rhf = scf.RHF(self.mol)
        rhf.kernel()
        c = rhf.mo_coeff
        hcore = rhf.get_hcore()
        hcore_mo = c.T @ hcore @ c
        self._num = hcore.shape[0]
        two_elec = self.mol.intor('int2e')
        # TWO_ELEC_MO = np.einsum('pi,qj,rk,sl,ijkl->pqrs', C, C, C, C, mol.intor('int2e'))
        two_elec_mo = ao2mo.kernel(self.mol, c, two_elec, compact=False)
        two_elec_mo = two_elec_mo.reshape((self._num, self._num, self._num, self._num))

        self._shots = None
        self.hamiltonian_pauli = self.hamiltonian(hcore_mo, two_elec_mo)

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

    def _circuit(self, qc, uccsd_ansatz):
        for qubit in range(self._num-1):
            qc.x(qubit)
            qc.x(qubit+self._num)

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

    def _measure(self, qc, shots):
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

            if all([p.symbol == 'I' for p in p_string.values()]):
                energy += values.real
                # print(f'Expectation: {values.real:18.15f} {p_string}')
                continue

            backend = AerProvider().get_backend('qasm_simulator')
            result = backend.run(qc_2, shots=shots).result().get_counts()

            counts = 0
            for key, value in result.items():
                counts += (-1)**sum([int(k) for k in key]) * value

            expectation = counts/shots * values.real

            # print(f'Expectation: {expectation:18.15f} {p_string}')
            energy += expectation
        return energy

    def _batch(self, coeff):
        uccsd_ansatz = self.uccsd_ansatz(coeff)
        qc = QuantumCircuit(2*self._num, 2*self._num)
        self._circuit(qc, uccsd_ansatz)
        energy = self._measure(qc, self._shots)
        print('coeff: ', end='')
        print([f'{val:9.6f}' for val in coeff], end='')
        print(f' Energy: {energy:18.15f}')
        return energy

    def run(self, shots=10000):
        self._shots = shots
        n = self._num
        coeff = [1e-5] * ((2 * (n//2) **2) + ((n//2)**2 * n * (n - 1) // 2))
        print(coeff)
        optimized_energy = opt.minimize(self._batch, coeff, method='Powell')

        nuclear_repulsion = self.mol.energy_nuc()
        total_energy = optimized_energy.fun + nuclear_repulsion
        print(f'Optimized Electronic Energy: {total_energy:18.15f}')
        return total_energy, optimized_energy.x
