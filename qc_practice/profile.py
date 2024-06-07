import json
import numpy as np
from qiskit import QuantumCircuit

class Profile:
    def __init__(self):
        self.state: int = 0
        self.num_orb: int = 0
        self.num_elec: int = 0
        self.energy_elec: float = 0.00
        self.energy_nucl: float = 0.00
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

    def copy(self):
        profile = Profile()
        profile.state = self.state
        profile.spin = self.spin
        profile.num_orb = self.num_orb
        profile.num_elec = self.num_elec
        profile.coeff = self.coeff
        profile.energy_elec = self.energy_elec
        profile.energy_nucl = self.energy_nucl
        profile.circuit = self.circuit

        return profile

    def save(self, filename):
        self.coeff = self.coeff.tolist()
        self.circuit = 'QuantumCircuit' # type: ignore
        with open(filename + '.json', 'w', encoding='utf-8') as f:
            json.dump(self.show(), f, indent=4)

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
    def __init__(self):
        self.profiles = []

    def __iter__(self):
        return iter(self.profiles)

    def __getitem__(self, key):
        return self.profiles[key]

    def __len__(self):
        return len(self.profiles)

    def add(self, profile, nroots):
        for _ in range(nroots):
            self.profiles.append(profile.copy())

    def show(self):
        return [profile.show() for profile in self.profiles]

    def energy_total(self):
        return [profile.energy_total() for profile in self.profiles]

    def save(self, filename):
        for i, profile in enumerate(self.profiles):
            profile.coeff = profile.coeff.tolist()
            profile.circuit = f'QuantumCircuit_{i}'
        with open(filename + '.json', 'w', encoding='utf-8') as f:
            json.dump(self.show(), f, indent=4)

    def __repr__(self):
        return f'{self.profiles}'
