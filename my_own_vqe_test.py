 #%%
from pyscf import gto
from qc_practice import VQE

mol = gto.M(atom = 'Li 0 0 0; H 0 0 1.596', basis = 'sto-3g')
vqe = VQE(mol)
vqe.run()

#%%
from pyscf import gto
from qc_practice import VQE
from qc_practice.ansatz import UCCSD

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.75', basis = 'sto-3g')
vqe = VQE(mol, UCCSD())
vqe.run()

#%%

from pyscf import gto
from qc_practice import SSVQE
from qc_practice.ansatz import UCCSD

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.74', basis = '4-31g')
ssvqe = SSVQE(mol, ansatz=UCCSD())
ssvqe.run()

#%%

from pyscf import gto
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider
from qc_practice import VQE


mol = gto.M(atom = 'H 0 0 0; H 0 0 0.74', basis = 'sto-3g')

vqe = VQE(mol)
e, coeff = vqe.run()
print(e, coeff)

for i in range(5):
    qc = QuantumCircuit(4, 4)
    vqe._initialize_circuit(qc)

    if i == 0:
        print('\nSinglet configuration 0110')
        qc.swap(0, 1)
    elif i == 1:
        print('\nSinglet configuration 1001')
        qc.swap(2, 3)
    elif i == 2:
        print('\nTriplet configuration 0011')
        qc.swap(0, 3)
    elif i == 3:
        print('\nTriplet configuration 1100')
        qc.swap(2, 1)
    else:
        print('\nSinglet configuration 1010')
        qc.swap(0, 1)
        qc.swap(2, 3)

    vqe._circuit(qc, vqe.uccsd_ansatz(coeff))
    energy = vqe._measure(qc) + mol.energy_nuc()

    print(energy)
    print('Spin measurement')

    qc.barrier()
    for k in range(4):
        qc.measure(k, k)

    # qc.draw('mpl')

    shots = 100000
    backend = AerProvider().get_backend('qasm_simulator')
    result = backend.run(qc, shots=shots).result().get_counts()
    print(result)

    spin = 0
    for key, value in result.items():
        for i, orb in enumerate(reversed([*key])):
            if i < 2:
                spin += value * int(orb)
            else:
                spin -= value * int(orb)

    print(f'multiplicity: {abs(spin/shots/2)}')
