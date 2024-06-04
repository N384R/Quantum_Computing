#%%
import json
import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto
from qc_practice import SSVQE
from qc_practice.ansatz import SpinFlipUCCSD

bond_lengths = np.arange(0.5, 2.5, 0.1)
potential_energy_surface = {}

for i in bond_lengths:
    mol = gto.M(atom = f'H 0 0 0; H 0 0 {i}', basis = 'sto-3g')
    ssvqe = SSVQE(mol)
    ssvqe.ansatz = SpinFlipUCCSD()
    ssvqe.weights = [1, 0.2, 0.2, 0.2, 0.05, 0.01]
    ssvqe.verbose = 0
    potential_energy_surface[f'{i:.02f}'] = ssvqe.run()
    print(f'\nBond Length: {i:5.2f}')
    print('Energy:')
    print(potential_energy_surface[f'{i:.02f}'])

plt.plot(list(potential_energy_surface.keys()), list(potential_energy_surface.values()))
plt.xlabel('Bond Length')
plt.ylabel('Energy')

plt.savefig('figures/H2_sfuccsd_sto3g_PES.png')
with open('data/H2_sfuccsd_sto3g_PES.json', 'w', encoding='utf-8') as f:
    json.dump(potential_energy_surface, f, indent=4)
