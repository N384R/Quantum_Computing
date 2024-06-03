class Profile:
    def __init__(self):
        self.state = None
        self.num_orb = None
        self.num_elec = None
        self.coeff = None
        self.energy_elec = None
        self.energy_nucl = None
        self.circuit = None
        self.spin = None

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
        for i in range(nroots):
            self.profiles.append(profile.copy())

    def show(self):
        return [profile.show() for profile in self.profiles]

    def __repr__(self):
        return f'{self.profiles}'
