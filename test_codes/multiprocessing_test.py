#%%
import time
from multiprocessing import Pool
import numpy as np
from pyscf import gto
from qc_practice import VQE
from qc_practice.ansatz import UCCSD

def pes_calc(bond_length):
    bl = f'{bond_length:.02f}'
    mol = gto.M(atom = f'H 0 0 0; H 0 0 {bond_length}', basis = 'sto-3g')
    vqe = VQE(mol)
    vqe.ansatz = UCCSD()
    vqe.verbose = 0
    energy = vqe.run().energy_total()
    return (bl, energy)

def pes_mult(bond_lengths):
    pes = {}
    with Pool(processes=4) as pool:
        for result in pool.imap(pes_calc, bond_lengths):
            bl, energy = result
            print(f'\nBond Length: {bl}')
            print(f'Energy: {energy}')
            pes[bl] = energy
    return pes

def non_parallel(bond_lengths):
    pes = {}
    for bl in bond_lengths:
        mol = gto.M(atom = f'H 0 0 0; H 0 0 {bl}', basis = 'sto-3g')
        vqe = VQE(mol)
        vqe.ansatz = UCCSD()
        vqe.verbose = 0
        energy = vqe.run().energy_total()
        print(f'\nBond Length: {bl:.02f}')
        print(f'Energy: {energy}')
        pes[f'{bl:.02f}'] = energy
    return pes

bond_lengths = np.arange(0.4, 1.0, 0.1)
start = int(time.time())
my_pes = pes_mult(bond_lengths)
end = int(time.time())
print(f"Time taken: {end-start} seconds for multiple process")

# start = int(time.time())
# my_pes = non_parallel(bond_lengths)
# end = int(time.time())
# print(f"Time taken: {end-start} seconds for single processes")
