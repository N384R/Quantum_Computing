#%%
import json
from multiprocessing import Pool
import numpy as np
import matplotlib.pyplot as plt
from pyscf import gto
from qc_practice import SSVQE
from qc_practice.ansatz import SpinFlipUCCSD

save_dir = '4-31g_test'

def ssvqe_calc(bond_len):
    mol = gto.M(atom = f'H 0 0 0; H 0 0 {bond_len}', basis = '4-31g')
    ssvqe = SSVQE(mol)
    ssvqe.ansatz = SpinFlipUCCSD()
    ssvqe.weights = [1, 0.1, 0.1, 0.1, 0.01, 0.001]
    ssvqe.verbose = 0
    ssvqe.run()
    return (ssvqe.profile, f'{bond_len:.02f}')

def pes_mult(bond_lens):
    surface = {}
    with Pool(processes=4) as pool:
        for result in pool.imap(ssvqe_calc, bond_lens):
            profile, bond_len = result
            bond_len = f'{bond_len}'
            print(f'\nBond Length: {bond_len}')
            energies = profile.energy_total()
            print(f'Energy: {energies}')
            surface |= {bond_len: energies}
            profile.save(f'{save_dir}/H2_sfuccsd_{bond_len}')
    return surface

bond_lengths = np.arange(0.50, 2.50, 0.1)
pes = pes_mult(bond_lengths)

x = np.array(list(pes.keys()), dtype=float)
data_np = np.array(list(pes.values()))
plt.plot(x, data_np)
plt.xlabel('Bond Length')
plt.ylabel('Energy')

plt.savefig(f'{save_dir}/H2_sfuccsd_PES.png')
with open(f'{save_dir}/H2_sfuccsd_PES.json', 'w', encoding='utf-8') as f:
    json.dump(pes, f, indent=4)
