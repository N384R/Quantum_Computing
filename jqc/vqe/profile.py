import json
from copy import deepcopy
import numpy as np
from qiskit import QuantumCircuit
from jqc.measure.hamiltonian import pyscf_luncher

class Profile:
    '''
    Class for storing the profile information of a molecule.
    '''
    def __init__(self, mol):
        self.info = pyscf_luncher(mol)

        self._state: int = 0
        self._spin: float = 0.00
        self._coeff: float = 0.00
        self._circuit: QuantumCircuit = QuantumCircuit()

    @property
    def num_orb(self) -> int:
        'Returns the number of orbitals.'
        return self.info['num_orb']

    @property
    def num_elec(self) -> int:
        'Returns the number of electrons.'
        return self.info['num_elec']

    @property
    def energy_nucl(self) -> float:
        'Returns the nuclear energy.'
        return self.info['energy_nuc']

    @property
    def energy_elec(self) -> float:
        'Returns the electronic energy.'
        return self.info['energy_elec']

    @energy_elec.setter
    def energy_elec(self, energy: float):
        self.info['energy_elec'] = energy

    @property
    def energy_total(self) -> float:
        'Returns the total energy of the system.'
        return self.energy_elec + self.energy_nucl

    @property
    def coeff(self) -> float:
        'Returns the coefficient of the state.'
        return self._coeff

    @coeff.setter
    def coeff(self, coeff: float):
        self._coeff = coeff

    @property
    def state(self) -> int:
        'Returns the state number.'
        return self._state

    @state.setter
    def state(self, state: int):
        self._state = state

    @property
    def spin(self) -> float:
        'Returns the spin of the state.'
        return self._spin

    @spin.setter
    def spin(self, spin: float):
        self._spin = spin

    @property
    def circuit(self) -> QuantumCircuit:
        'Returns the quantum circuit'
        return self._circuit

    @circuit.setter
    def circuit(self, circuit: QuantumCircuit):
        self._circuit = circuit

    @property
    def qm(self) -> dict:
        'Returns the quantum mechanical information.'
        return {
            'hcore': self.info['hcore_mo'],
            'two_elec': self.info['two_elec_mo'],
            'dipole': self.info['dipole_mo']
        }

    def __getitem__(self, key):
        return getattr(self, key)

    def show(self):
        'Returns the profile information.'
        return {
            'state': self.state,
            'spin': f'{self.spin:.02f}',
            'num_orb': self.num_orb,
            'num_elec': self.num_elec,
            'coeff': self.coeff,
            'energy_total': self.energy_total,
            'circuit': self.circuit,
        }

    def save(self, filename):
        'Saves the profile information to a JSON file.'
        def convert(o):
            if isinstance(o, np.ndarray):
                return o.tolist()
            if isinstance(o, QuantumCircuit):
                return 'QuantumCircuit'
            return o
        with open(filename + '.json', 'w', encoding='utf-8') as f:
            json.dump(self.show(), f, indent=4, default=convert)

    def __repr__(self):
        if round(self.spin) == 0:
            mult = 'Singlet'
        elif round(self.spin) == 2:
            mult = 'Triplet'
        else:
            mult = 'Unknown'
        state_number = self.state if self.state else 0
        return f'{mult}_{state_number}'

class Profiles:
    '''
    Class for storing the profile information of a molecule.
    '''
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
        'Update the profile with the new state.'
        self.profiles[profile1.state] = deepcopy(profile1)

    def show(self):
        'Returns the profile information.'
        return [profile.show() for profile in self.profiles]

    def energy_total(self):
        'Returns the total energy of the states.'
        return [profile.energy_total() for profile in self.profiles]

    def __repr__(self):
        return f'{self.profiles}'
