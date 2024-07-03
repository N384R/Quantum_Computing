import json
from copy import deepcopy
import numpy as np
from qiskit import QuantumCircuit
from qc_practice.measure.hamiltonian import pyscf_luncher

class Profile:
    def __init__(self, mol):
        info = pyscf_luncher(mol)

        self.num_orb = info['num_orb']
        self.num_elec = info['num_elec']
        self.energy_nucl = info['energy_nuc']
        self.energy_elec = info['energy_elec']
        self.qm = {'hcore': info['hcore_mo'], 'two_elec': info['two_elec_mo']}

        self.state: int = 0
        self.spin: float = 0.00
        self.coeff: np.ndarray = np.array([])
        self.circuit: QuantumCircuit = QuantumCircuit()

    def energy_total(self):
        return self.energy_elec + self.energy_nucl

    def show(self):
        return {
            'state': self.state,
            'spin': f'{self.spin:.02f}',
            'num_orb': self.num_orb,
            'num_elec': self.num_elec,
            'coeff': self.coeff,
            'energy_elec': self.energy_elec,
            'energy_nucl': self.energy_nucl,
            'circuit': self.circuit,
        }

    def save(self, filename):
        def convert(o):
            if isinstance(o, np.ndarray):
                return o.tolist()
            if isinstance(o, QuantumCircuit):
                return 'QuantumCircuit'
            return o
        with open(filename + '.json', 'w', encoding='utf-8') as f:
            json.dump(self.show(), f, indent=4, default=convert)

    #Need to me modified to print currect state
    #Need to me modified to print currect state
    #Need to me modified to print currect state
    def __repr__(self):
        if round(self.spin) == 0:
            mult = 'Singlet'
        elif round(self.spin) == 1:
            mult = 'Triplet'
        else:
            mult = 'Unknown'
        state_number = self.state if self.state else 0
        return f'{mult}_{state_number}'

class Profiles:
    def __init__(self, p, nstates=1):
        self.profiles = self._generate_profiles(p, nstates)

    def __iter__(self):
        return iter(self.profiles)

    def __getitem__(self, key):
        return self.profiles[key]

    def __len__(self):
        return len(self.profiles)

    @staticmethod
    def _generate_profiles(p, nstates):
        profiles = []
        for _ in range(nstates):
            profiles.append(deepcopy(p))
        return profiles

    def update(self, profile1):
        self.profiles[profile1.state] = deepcopy(profile1)

    def show(self):
        return [profile.show() for profile in self.profiles]

    def energy_total(self):
        return [profile.energy_total() for profile in self.profiles]

    def __repr__(self):
        return f'{self.profiles}'
