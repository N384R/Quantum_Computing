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

print("nuclear_repulsion_energy", hamiltonian.nuclear_repulsion_energy)n

print("reference_energy", problem.reference_energy)

print("num_particles:", problem.num_particles)

print("spatial_orbitals:", problem.num_spatial_orbitals)

