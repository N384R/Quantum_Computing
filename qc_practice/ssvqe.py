from math import ceil
from pyscf import scf, ao2mo
import scipy.optimize as opt
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import JordanWignerMapper
from qc_practice import VQE

class SSVQE(VQE):
    def __init__(self, mol, nroots = 2, weights = None):
        self.mol = mol
        self.verbose = 1
        self.nroots = nroots * 2
        self.state_energies = []
        self._current_state = None
        self._trial = 0

        if weights is None:
            self.weights = [1]
            for i in range(1, self.nroots):
                self.weights.append(self.weights[i-1] / 4)
        else:
            self.weights = weights

    def _batch(self, coeff):
        uccsd_ansatz = self.uccsd_ansatz(coeff)
        qc = QuantumCircuit(2*self._num, 2*self._num)
        self._initialize(qc)
        N = self._num
        CN = ceil(N/3)
        if self._current_state != 0:
            if self._current_state % 3 == 1:
                qc.swap(N//2-1, N + CN)
            elif self._current_state % 3 == 2:
                qc.swap(N//2-1, N//2 - 1 + CN)
            else:
                qc.swap(N//2-1, N//2 - 1 + CN)
                qc.swap(N, N + CN)

        self._circuit(qc, uccsd_ansatz)
        energy = self._measure(qc)

        if self._current_state % 3 == 0:
            self._talk('coeff: ', end='', verb=2)
            self._talk([f'{val:9.6f}' for val in coeff], verb=2)

        return energy

    def _ssvqe_batch(self, coeff):
        self._trial += 1
        self.state_energies = []
        self._talk(f'\nTrial: {self._trial}')
        for i in range(self.nroots):
            self._iterations = 0
            self._current_state = i
            state_energy = self._batch(coeff)
            self._talk(f'State_{i} Energy: {state_energy + self.nuclear_repulsion:18.15f}')
            self.state_energies.append(state_energy)
        weighted_energy = sum(self.weights[i] * self.state_energies[i] for i in range(self.nroots))

        return weighted_energy

    def run(self, shots=10000):
        super().__init__(self.mol, verbose=self.verbose)
        self._talk('Running SSVQE')
        n = self._num
        self._shots = shots
        count = (2 * (n//2) ** 2) + 2 * (n//2 * (n//2 - 1) // 2) ** 2 + (n//2) ** 4
        coeff = [1e-5] * count
        opt.minimize(self._ssvqe_batch, coeff, method='Powell')
        self._talk('\nFinal Excited State Energies:')
        for i in range(self.nroots):
            self.state_energies[i] = self.state_energies[i] + self.nuclear_repulsion
            self._talk(f'State_{i} Energy: {self.state_energies[i]:18.15f}', end='\n' if i == 0 else '  ')
            if i !=0:
                energy_diff = self.state_energies[i] - self.state_energies[0]
                self._talk(f'Excitation energy: {energy_diff:18.15f}')

        return self.state_energies
