from itertools import combinations
import scipy.optimize as opt
from qc_practice import VQE
from qc_practice.profile import Profiles
from qc_practice.ansatz.uccsd import UCCSD

class SSVQE(VQE):
    def __init__(self, mol, ansatz=UCCSD(), active_space = None, weights = None, koopmans=False, verbose=1):
        self.mol = mol
        self._profiles = Profiles()

        self.verbose = verbose
        self._transition = None
        self.koopmans = koopmans
        self._trial = 0

        if active_space is None:
            self.active_space = [1, 1]
        else:
            self.active_space = active_space
        k = self.active_space
        self._nspace = k[0] * k[1] * 4
        if not koopmans:
            self._nspace += k[0] * (2 * k[0] - 1) * k[1] * (2 * k[1] - 1)

        if weights is None:
            self.weights = [1]
            for i in range(1, self._nspace+1):
                self.weights.append(self.weights[i-1] / 2)
        else:
            self.weights = weights

        super().__init__(mol)
        if sum([*self.active_space]) > self.profile.num_orb:
            raise ValueError('Active space is larger than number of orbitals')

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

    def _swap_init(self, qc, state):
        state = next(self._transition)
        qc.swap(state[0], state[1])
        print(f'Swapping {state[0]} and {state[1]}')

    def _initialize_circuit(self, qc):
        super()._initialize_circuit(qc)
        if self.profile.state != 0:
            self._swap_init(qc, self.profile.state)

    def _ssvqe_batch(self, coeff):
        self._trial += 1
        self._talk(f'\nIteration: {self._trial}')
        self._transition = self._electron_excitation(self.active_space)
        for i in range(self._nspace+1):
            self.profile.state = i
            self._profiles[i].state = i
            self.verbose, verbose = 0, self.verbose
            self._profiles[i].energy_elec = self._batch(coeff)
            self._profiles[i].circuit = self.profile.circuit
            self.verbose = verbose
            self._talk(f'State_{i} Energy: {self._profiles[i].energy_total():18.15f}')
        if self._trial == 1:
            indices = sorted(range(len(self._profiles)), key = lambda x: self._profiles[x].energy_elec)
            self.weights = [self.weights[i] for i in indices]
        weighted_energy = sum(self.weights[i] * self._profiles[i].energy_elec for i in range(self._nspace+1))

        return weighted_energy

    def run(self, shots=10000):
        self._talk('\nStarting SSVQE Calculation')
        self._profiles.add(self.profile, self._nspace+1)
        n = self.profile.num_orb
        self._shots = shots
        coeff = [1e-5] * ((2 * (n//2) **2) + 2 * (n//2 * (n//2 - 1) // 2)**2 + (n//2)**4)
        optimized = opt.minimize(self._ssvqe_batch, coeff, method='Powell')
        self._talk('\n!!Successfully Converged!!')
        self.profile.coeff = optimized.x
        self.profile = self._profiles_update()
        self._talk('\nFinal Excited State Energies:')
        for i in range(self._nspace+1):
            self._talk(f'State_{i} Energy: {self._profiles[i].energy_total():18.15f}')

    def _profiles_update(self):
        for i in range(self._nspace+1):
            self._profiles[i].coeff = self.profile.coeff
            self._profiles[i].spin = self._measure_spin(self._profiles[i].circuit)
        return self._profiles
