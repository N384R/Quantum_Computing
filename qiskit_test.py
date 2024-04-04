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

result = solver.solve(problem)
print(result)

#%%

from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_algorithms import NumPyMinimumEigensolver

from qiskit_algorithms import VQE
from qiskit_algorithms.optimizers import SLSQP
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

es_problem = driver.run()
mapper = JordanWignerMapper()
numpy_solver = NumPyMinimumEigensolver()

ansatz = UCCSD(
    es_problem.num_spatial_orbitals,
    es_problem.num_particles,
    mapper,
    initial_state=HartreeFock(
        es_problem.num_spatial_orbitals, 
        es_problem.num_particles,
        mapper
    ),
)

vqe_solver = VQE(Estimator(), ansatz, SLSQP())
vqe_solver.initial_point = [0.0] * ansatz.num_parameters

calc = GroundStateEigensolver(mapper, vqe_solver)
res = calc.solve(es_problem)
print(res)

calc = GroundStateEigensolver(mapper, numpy_solver)
res = calc.solve(es_problem)
print(res)