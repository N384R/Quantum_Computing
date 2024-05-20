#%%
from pyscf import gto, scf, ao2mo
from qiskit import QuantumCircuit
from qiskit_aer import AerProvider

from qc_practice import JordanWignerMapper

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.735', basis = 'sto-3g')
rhf = scf.RHF(mol)
EN = rhf.kernel()
C = rhf.mo_coeff
HCORE = rhf.get_hcore()
HCORE_MO = C.T @ HCORE @ C
S = rhf.get_ovlp()
OCC = rhf.get_occ()
NUM = HCORE.shape[0]
TWO_ELEC = mol.intor('int2e')
# TWO_ELEC_MO = np.einsum('pi,qj,rk,sl,ijkl->pqrs', C, C, C, C, mol.intor('int2e'))
TWO_ELEC_MO = ao2mo.kernel(mol, C, TWO_ELEC, compact=False)
TWO_ELEC_MO = TWO_ELEC_MO.reshape((NUM, NUM, NUM, NUM))
second_q = ''
for i in range(NUM):
    coeff = HCORE_MO[i, i]
    SIGN = '+' if coeff > 0 else ''
    second_q += f'{SIGN} {coeff:.16f} {i}^ {i}' + '\n'
    second_q += f'{SIGN} {coeff:.16f} {i+NUM}^ {i+NUM}' + '\n'

two_elec_dict = {}
for i in range(NUM):
    for j in range(NUM):
        for k in range(NUM):
            for l in range(NUM):
                coeff = TWO_ELEC_MO[i, j, k, l]/2
                if coeff < 1e-10:
                    continue
                SIGN = '+' if coeff > 0 else ''
                second_q += f'{SIGN} {coeff:.16f} {i}^ {k}^ {l} {j}' + '\n'
                second_q += f'{SIGN} {coeff:.16f} {i}^ {k+NUM}^ {l+NUM} {j}' + '\n'
                second_q += f'{SIGN} {coeff:.16f} {i+NUM}^ {k}^ {l} {j+NUM}' + '\n'
                second_q += f'{SIGN} {coeff:.16f} {i+NUM}^ {k+NUM}^ {l+NUM} {j+NUM}' + '\n'

hamiltonian_pauli = JordanWignerMapper(second_q)
print(hamiltonian_pauli)

coeff = 1e-5
uccsd_ansatz = ''
for i in range(NUM//2):
    for j in range(NUM//2, NUM):
        uccsd_ansatz += f'+ {coeff:f} {i}^ {j} ' + '\n'
        uccsd_ansatz += f'- {coeff:f} {j}^ {i} ' + '\n'
        uccsd_ansatz += f'+ {coeff:f} {i+NUM}^ {j+NUM} ' + '\n'
        uccsd_ansatz += f'- {coeff:f} {j+NUM}^ {i+NUM} ' + '\n'

for i in range(NUM//2):
    for j in range(NUM//2):
        for k in range(NUM//2, NUM):
            for l in range(NUM//2, NUM):
                uccsd_ansatz += f'+ {coeff:f} {i}^ {j+NUM}^ {k} {l+NUM} ' + '\n'
                uccsd_ansatz += f'- {coeff:f} {k}^ {l+NUM}^ {i} {j+NUM} ' + '\n'

uccsd_pauli = JordanWignerMapper(uccsd_ansatz)
print(uccsd_pauli)

qc = QuantumCircuit(2*NUM, 2*NUM)
for qubit in range(NUM-1):
    qc.x(qubit)
    qc.x(qubit+NUM)

for p_string, values in uccsd_pauli.items():

    chk = []
    for idx, p in p_string.items():
        if p.symbol == 'X':
            qc.h(idx)

        elif p.symbol == 'Y':
            qc.s(idx)
            qc.h(idx)

        if p.symbol != 'I':
            chk += [idx]

    for i in chk:
        if i != max(chk):
            qc.cx(i, i+1)

    qc.rz(values.imag, max(chk))

    for i in reversed(chk):
        if i != max(chk):
            qc.cx(i, i+1)

    for idx, p in p_string.items():
        if p.symbol == 'X':
            qc.h(idx)

        elif p.symbol == 'Y':
            qc.h(idx)
            qc.sdg(idx)

qc.barrier()

def measure(qc, hamiltonian_pauli, shots=10000):
    total_energy = 0.
    for p_string, values in hamiltonian_pauli.items():
        for idx, p in p_string.items():
            if p.symbol == 'X':
                qc.h(idx)

            elif p.symbol == 'Y':
                qc.sdg(idx)
                qc.h(idx)

        qc.barrier()

        for idx, p in p_string.items():
            if p.symbol == 'I':
                continue

            qc.measure(idx, idx)
        
        qc.barrier()

        if all([p.symbol == 'I' for p in p_string.values()]):
            total_energy += values.real
            print(f'Expectation: {values.real:18.15f} {p_string}')
            continue

        backend = AerProvider().get_backend('qasm_simulator')
        result = backend.run(qc, shots=shots).result().get_counts()

        counts = 0
        for key, value in result.items():
            counts += (-1)**sum([int(k) for k in key]) * value

        expectation = counts/shots * values.real

        print(f'Expectation: {expectation:18.15f} {p_string}')
        total_energy += expectation
    return total_energy

total_energy = measure(qc, hamiltonian_pauli)
print(total_energy)


qc.draw('mpl')
