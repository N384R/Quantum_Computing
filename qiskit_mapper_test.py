from qiskit_algorithms import NumPyMinimumEigensolver
from qiskit_nature.second_q.algorithms import GroundStateEigensolver
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit

driver = PySCFDriver(
    atom = "H 0 0 0; H 0 0 0.735",
    basis = 'sto3g',
    charge = 0,
    spin = 0,
    unit = DistanceUnit.ANGSTROM,
)

problem = driver.run()
fermionic_op = problem.hamiltonian.second_q_op()
mapper = JordanWignerMapper()
qubit_jw_op = mapper.map(fermionic_op)
print(fermionic_op)
print(qubit_jw_op)

