from typing import cast
from itertools import combinations
import scipy.optimize as opt
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import SSVQE
from .profile import Profiles

class VQD(SSVQE):
    def __init__(self, mol, ansatz = None):
        super().__init__(mol, ansatz = ansatz)

    def _excite_batch(self, coeff):
        self._config['trial'] += 1
        self._talk(f"\nIteration: {self._config['trial']}")
        self._transition = self._electron_excitation(self.active_space)
        for i in range(self._config['nspace']+1):
            self._config['state'] = i
            self.verbose, verbose = 0, self.verbose
            self._p[i].energy_elec = self._batch(coeff)
            self._p[i].transition = self._transition
            self._p[i].circuit = self.profile.circuit
            self._p[i].state = i
            self.verbose = verbose
            self._talk(f'State_{i} Energy: {self._p[i].energy_total():18.15f}')
        weight_en = sum(self.weights[i] * self._p[i].energy_elec
                        for i in range(self._config['nspace']+1))

        return weight_en

    def _calculate_overlap(self):
        pass
