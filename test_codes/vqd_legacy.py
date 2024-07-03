from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import VQE
from ..qc_practice.profile import Profiles

class VQD(VQE):
    def __init__(self, mol, ansatz = None, simulator = None):
        super().__init__(mol, ansatz = ansatz, simulator = simulator)
        self._config['trial'] = 0
        self._config['nstates'] = 2
        self._config['verbose'] -= 1

    @property
    def nstates(self):
        'The number of states for the calculation.'
        return self._config['nstates']

    @nstates.setter
    def nstates(self, nstates):
        self._config['nstates'] = nstates

    def _excite_batch(self, coeff):
        beta = 3.0
        self.profile.coeff = coeff
        self.profile.energy_elec = self._batch(coeff)
        energy = self.profile.energy_elec
        for i in range(self.profile.state):
            state1 = self.profile
            state2 = self._config['profiles'][i]
            overlap_sq = self._measure_overlap(state1, state2)
            energy += beta * overlap_sq
        return energy

    def run(self, shots=10000):
        self._introduce()
        self._init_setup()
        self._config['shots'] = shots
        self._config['profiles'] = Profiles(self.profile, self._config['nstates'])

        self._talk('\nGround state calculation', verb=1)
        self._talk("\nState: 0", verb=1)
        coeff = self.ansatz.generate_coeff(self.profile)
        optimized = self.ansatz.call_optimizer(self._batch, coeff, self.optimizer)

        self.profile.coeff = optimized.x
        self.profile.energy_elec = optimized.fun
        self._config['profiles'].update(self.profile)
        self._result()

        self._talk('\nExcited state calculation', verb=1)
        for state in range(1, self._config['nstates']):
            self.profile.state = state
            self._config['iteration'] = 0
            self._talk(f"\nState: {state}", verb=1)
            coeff = self.ansatz.generate_coeff(self.profile)
            optimized = self.ansatz.call_optimizer(self._excite_batch, coeff, self.optimizer)

            self.profile.coeff = optimized.x
            self.profile.energy_elec = optimized.fun
            self._config['profiles'].update(self.profile)
            self._result()

        self._talk('\n!!Calculation Completed!!\n', verb=1)
        self._talk('Final Excited State Energies:', verb=1)
        for i in range(self._config['nstates']):
            self._talk(f"State_{i} Energy: {self._config['profiles'][i].energy_total():18.15f}"
                        , verb=1)
        return self.profile


    def _measure_overlap(self, state1, state2):
        return self._swap_test(state1, state2)

    def _swap_test(self, state1, state2):
        no = self.profile.num_orb
        ancila = QuantumCircuit(1, 1)
        qc1 = QuantumCircuit(no*2)
        self._initialize_circuit(qc1)
        self.ansatz.ansatz(qc1, self.profile, state1.coeff)         # need to modify by using: state1.circuit

        qc2 = QuantumCircuit(no*2)
        self._initialize_circuit(qc2)
        self.ansatz.ansatz(qc2, self.profile, state2.coeff)

        qc3 = qc1.tensor(qc2).tensor(ancila)
        qc3.h(0)
        for i in range(1, no*2+1):
            qc3.cswap(0, i, i+no*2)
        qc3.h(0)
        qc3.measure(0, 0)

        backend = AerProvider().get_backend('qasm_simulator')
        result = backend.run(qc3, shots=self._config['shots']).result().get_counts()
        overlap_squared = abs(result.get('0') / self._config['shots'] * 2 - 1)
        return overlap_squared

    def _introduce(self):
        self._talk(f'\nStarting {self.__class__.__name__} Calculation\n', verb=1)
        self._talk(f'Ansatz: {self.ansatz.__class__.__name__}', verb=1)
        self._talk(f'Optimizer: {self.optimizer}', verb=1)
        self._talk(f'Number of states: {self.nstates}', verb=1)

    def _result(self):
        self._talk('\n\n!!Successfully Converged!!', verb=1)
        self._talk(f"Optimized Electronic Energy: {self.profile.energy_elec:18.15f}", verb=1)
        self._talk(f"Optimized Total Energy     : {self.profile.energy_total():18.15f}", verb=1)
