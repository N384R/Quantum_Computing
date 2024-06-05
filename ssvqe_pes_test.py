#%%
import json
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto
from qc_practice import SSVQE
from qc_practice.ansatz import SpinFlipUCCSD

def ssvqe_calc(bond_len):
    mol = gto.M(atom = f'H 0 0 0; H 0 0 {bond_len}', basis = 'sto-3g')
    ssvqe = SSVQE(mol)
    ssvqe.ansatz = SpinFlipUCCSD()
    ssvqe.weights = [1, 0.5, 0.5, 0.5, 0.1, 0.01]
    ssvqe.verbose = 0
    ssvqe.run()
    return (ssvqe.profile, f'{bond_len}')

def pes_mult(bond_lens):
    surface = {}
    with Pool(processes=4) as pool:
        for result in pool.imap(ssvqe_calc, bond_lens):
            profile, bond_len = result
            print(f'\nBond Length: {bond_len}')
            energies = profile.energy_total()
            print(f'Energy: {energies}')
            surface |= {bond_len: energies}
            profile.save(f'pes_test/H2_sto3g_{bond_len}')
    return surface

bond_lengths = np.arange(0.5, 2.5, 0.1)
pes = pes_mult(bond_lengths)

plt.plot(list(pes.keys()), list(pes.values()))
plt.xlabel('Bond Length')
plt.ylabel('Energy')

plt.savefig('pes_test/H2_sfuccsd_sto3g_PES.png')
with open('pes_test/H2_sfuccsd_sto3g_PES.json', 'w', encoding='utf-8') as f:
    json.dump(pes, f, indent=4)
