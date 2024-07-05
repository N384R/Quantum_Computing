#%%
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver

driver = PySCFDriver(
    atom="H 0 0 0; H 0 0 0.735",
    basis="sto3g",
    charge=0,
    spin=0,
    unit=DistanceUnit.ANGSTROM,
)

problem = driver.run()

hamiltonian = problem.hamiltonian
coefficients = hamiltonian.electronic_integrals
print(coefficients.alpha)

second_q_op = hamiltonian.second_q_op()
print(second_q_op)

print("nuclear_repulsion_energy", hamiltonian.nuclear_repulsion_energy)
print("reference_energy", problem.reference_energy)
print("num_particles:", problem.num_particles)
print("spatial_orbitals:", problem.num_spatial_orbitals)

#%%

from qiskit_algorithms import NumPyMinimumEigensolver
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit

driver = PySCFDriver(
    atom="H 0 0 0; H 0 0 0.735",
    basis="sto3g",
    charge=0,
    spin=0,
    unit=DistanceUnit.ANGSTROM,
)

problem = driver.run()

solver = GroundStateEigensolver(
    JordanWignerMapper(),
    NumPyMinimumEigensolver(),
)

result_singlet = solver.solve(problem)
print(result_singlet)

#%%

from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_algorithms import NumPyMinimumEigensolver, NumPyEigensolver

from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import POWELL
from qiskit.primitives import Estimator
from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD

from qiskit_nature.second_q.algorithms import GroundStateEigensolver

driver = PySCFDriver(
    atom = "H 0 0 0; H 0 0 0.735",
    basis = 'sto3g',
    charge = 0,
    spin = 0,
    unit = DistanceUnit.ANGSTROM,
)

problem_singlet = driver.run()
second_q_op = problem_singlet.hamiltonian.second_q_op()
mapper = JordanWignerMapper()
numpy_solver = NumPyMinimumEigensolver()
qubit_op = mapper.map(second_q_op)

ansatz = UCCSD(
    problem_singlet.num_spatial_orbitals,
    problem_singlet.num_particles,
    mapper,
    initial_state=HartreeFock(
        problem_singlet.num_spatial_orbitals,
        problem_singlet.num_particles,
        mapper
    ),
)

ansatz.decompose().decompose().draw('mpl')

vqe_solver = VQE(Estimator(), ansatz, POWELL())
vqe_solver.initial_point = [0.0] * ansatz.num_parameters

calc = GroundStateEigensolver(mapper, vqe_solver)
res = calc.solve(problem_singlet)
print(res)

# calc = GroundStateEigensolver(mapper, numpy_solver)
# res = calc.solve(es_problem)
# print(res)

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
    print('singlet', result_singlet.total_energies)
    potential_energy_surface[f'{i:.02f}'] = list(result_singlet.total_energies)

    # driver_triplet = PySCFDriver(
    #     atom = f"H 0 0 0; H 0 0 {i}",
    #     basis = 'sto3g',
    #     charge = 0,
    #     spin = 2,
    #     unit = DistanceUnit.ANGSTROM,
    #     method=MethodType.ROHF,
    # )

    # problem_triplet = driver_triplet.run()
    # mapper = JordanWignerMapper()
    # es_solver_triplet = NumPyEigensolver(k=3)
    # es_solver_triplet.filter_criterion = problem_triplet.get_default_filter_criterion()
    # solver = ExcitedStatesEigensolver(mapper, es_solver_triplet)

    # result_triplet = solver.solve(problem_triplet)
    # print('triplet', result_triplet.total_energies)
    # potential_energy_surface[i].extend(result_triplet.total_energies)

plt.plot(list(potential_energy_surface.keys()), list(potential_energy_surface.values()))
plt.xlabel('Bond Length')
plt.ylabel('Energy')
plt.savefig('figures/H2_exact_PES.png')

with open('data/H2_exact_PES.json', 'w', encoding='utf-8') as f:
    json.dump(potential_energy_surface, f, indent=4)
