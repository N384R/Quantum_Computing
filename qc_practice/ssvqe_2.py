from pyscf import scf, ao2mo
import scipy.optimize as opt
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import JordanWignerMapper
from qc_practice import VQE

class SSVQE(VQE):
    def __init__(self, mol, nstates = 2, weights = [1, 1]):
        super().__init__(mol)
        self.mol = mol
        self.nstates = nstates
        self.weights = weights
        self.energies = []
        self._current_state = None
        self._trial = 0

    def _initialize(self, qc):
        for qubit in range(self._num_elec//2):
            qc.x(qubit+1)
            qc.x(qubit+self._num)

    def _vqe_batch(self, coeff):
        self._iterations += 1
        print(f'State_{self._current_state} iteration: {self._iterations}', end='\r', flush=True)
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

    def _vqe_run(self, coeff):
        optimized_energy = opt.minimize(self._vqe_batch, coeff, method='Powell')
        return optimized_energy.fun

    def _ssvqe_batch(self, coeff):
        self._trial += 1
        state_energies = []
        print(f'Trial: {self._trial}')
        for i in range(self.nstates):
            self._iterations = 0
            self._current_state = i
            state_energy = self._vqe_run(coeff)
            print(f'State_{i} Energy: {state_energy + self.nuclear_repulsion:18.15f}')
            state_energies.append(state_energy)
        weighted_energy = sum(self.weights[i] * state_energies[i] for i in range(self.nstates))
        return weighted_energy

    def run(self, shots=10000):
        print('Running SSVQE\n')
        n = self._num
        self._shots = shots
        coeff = [1e-5] * ((2 * (n//2) **2) + 2 * (n//2 * (n//2 - 1) // 2)**2 + (n//2)**4)
        optimized_excited_states = opt.minimize(self._ssvqe_batch, coeff, method='Powell')
        value, contributions = optimized_excited_states.fun, optimized_excited_states.x
        print('Final Excited State Energies:')
        for i in range(self.nstates):
            state_energy = value / contributions[i] + self.nuclear_repulsion
            print(f'State_{i} Energy: {state_energy:18.15f}')
            self.energies.append(state_energy)
        return self.energies
