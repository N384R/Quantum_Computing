#%%

from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver, MethodType
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_algorithms import NumPyEigensolver
from qiskit_nature.second_q.algorithms import ExcitedStatesEigensolver

from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import POWELL
from qiskit.primitives import Estimator
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
import numpy as np
import json
import matplotlib.pyplot as plt
from pyscf import gto, scf
from pyscf import cc

exact = True

bond_lengths = np.arange(0.5, 2.5, 0.1)
potential_energy_surface = {}

for i in bond_lengths:
    driver_singlet = PySCFDriver(
        atom = f"H 0 0 0; H 0 0 {i}",
        basis = 'sto3g',
        charge = 0,
        spin = 0,
        unit = DistanceUnit.ANGSTROM,
        method=MethodType.ROHF,
    )

    problem_singlet = driver_singlet.run()
    mapper = JordanWignerMapper()
    es_solver_singlet = NumPyEigensolver(k=3)
    es_solver_singlet.filter_criterion = problem_singlet.get_default_filter_criterion()
    solver = ExcitedStatesEigensolver(mapper, es_solver_singlet)

    result_singlet = solver.solve(problem_singlet)
    # print('singlet', result_singlet.total_energies)
    potential_energy_surface[f'{i:.02f}'] = list(result_singlet.total_energies)

for key, value in potential_energy_surface.items():
    mol = gto.M(atom = f'H 0 0 0; H 0 0 {key}', basis = 'sto-3g', spin=2)
    mf = scf.UHF(mol)
    mf.kernel()
    ccsd = cc.CCSD(mf)
    e = ccsd.kernel()
    total_e = mf.e_tot + e[0] if exact else mf.e_tot
    potential_energy_surface[key].insert(1, total_e)
    print(f'\nBond Length: {key}')
    print('Energy:')
    print(value)

data_np = np.array(list(potential_energy_surface.values()))
x = np.array(list(potential_energy_surface.keys()), dtype=float)
plt.plot(x, data_np)
plt.xlabel('Bond Length')
plt.ylabel('Energy')
plt.savefig('figures/H2_exact_PES.png')

with open('data/H2_exact_PES.json', 'w', encoding='utf-8') as f:
    json.dump(potential_energy_surface, f, indent=4)
