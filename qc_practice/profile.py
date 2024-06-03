class Profile:
    def __init__(self):
        self.num_orb = None
        self.num_elec = None
        self.coeff = None
        self.energy_elec = None
        self.energy_nucl = None
        self.circuit = None

    def energy_total(self):
        return self.energy_elec + self.energy_nucl

    def show(self):
        return {
            'num_orb': self.num_orb,
            'num_elec': self.num_elec,
            'coeff': self.coeff,
            'energy_elec': self.energy_elec,
            'energy_nucl': self.energy_nucl,
            'circuit': self.circuit,
        }
