from multiprocessing import Pool
from qiskit import QuantumCircuit
from qc_practice.simulator import Simulator

def measure(qc: QuantumCircuit, h_pauli, simulator: Simulator, parallel: bool = False):
    'Measure the expectation value of a Hamiltonian'
    tasks = [(qc, p_string, values, simulator) for p_string, values in h_pauli.items()]

    if parallel:
        with Pool(4) as pool:
            energy = sum(list(pool.map(single_measure, tasks)))

    else:
        energy = sum(single_measure(task) for task in tasks)

    return energy

def single_measure(args: tuple[QuantumCircuit, dict, complex, Simulator]):
    'Measure the expectation value of a Pauli string'
    qc, p_string, values, simulator = args

    if all(p.symbol == 'I' for p in p_string.values()):
        return values.real

    qc_2 = qc.copy()
    probability = simulator.run_simulator(qc_2, p_string)

    expectation = probability * values.real
    return expectation
