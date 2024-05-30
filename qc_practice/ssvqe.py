from pyscf import scf, ao2mo
import scipy.optimize as opt
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import JordanWignerMapper
from qc_practice import VQE

class SSVQE(VQE):
    def __init__(self, mol, nstates = 2, weights = None):
        super().__init__(mol)
        self.nstates = nstates
        self.state_energies = []
        self._current_state = None
        self._trial = 0

        if weights is None:
            self.weights = [1]
            for i in range(1, self.nstates):
                self.weights.append(self.weights[i-1] / 2)
        else:
            self.weights = weights

    def _initialize(self, qc):
        for qubit in range(self._num_elec//2):
            qc.x(qubit+1)
            qc.x(qubit+self._num)

    def _vqe_batch(self, coeff):
        uccsd_ansatz = self.uccsd_ansatz(coeff)
        qc = QuantumCircuit(2*self._num, 2*self._num)
        self._initialize(qc)
        if self._current_state != 0:
            qc.swap(self._num//2-1, self._num//2 + self._current_state-1)
        self._circuit(qc, uccsd_ansatz)
        energy = self._measure(qc)

        if self.verbose:
            print('coeff: ', end='')
            print([f'{val:9.6f}' for val in coeff])
        return energy

    def _ssvqe_batch(self, coeff):
        self._trial += 1
        self.state_energies = []
        print(f'\nTrial: {self._trial}')
        for i in range(self.nstates):
            self._iterations = 0
            self._current_state = i
            state_energy = self._vqe_batch(coeff)
            print(f'State_{i} Energy: {state_energy + self.nuclear_repulsion:18.15f}')
            self.state_energies.append(state_energy)
        weighted_energy = sum(self.weights[i] * self.state_energies[i] for i in range(self.nstates))
        return weighted_energy

    def run(self, shots=10000):
        print('Running SSVQE')
        n = self._num
        self._shots = shots
        coeff = [1e-5] * ((2 * (n//2) **2) + 2 * (n//2 * (n//2 - 1) // 2)**2 + (n//2)**4)
        opt.minimize(self._ssvqe_batch, coeff, method='Powell')
        print('\nFinal Excited State Energies:')
        for i in range(self.nstates):
            state_energy = self.state_energies[i] + self.nuclear_repulsion
            print(f'State_{i} Energy: {state_energy:18.15f}')
        return self.state_energies
