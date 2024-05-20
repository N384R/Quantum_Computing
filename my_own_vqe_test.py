#%%
from pyscf import gto, scf, ao2mo
from qc_practice import JordanWignerMapper
from qiskit import QuantumCircuit

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
SECOND_Q = ''
for i in range(NUM):
    coeff = HCORE_MO[i, i]
    SIGN = '+' if coeff > 0 else ''
    SECOND_Q += f'{SIGN} {coeff:.16f} {i}^ {i}' + '\n'
    SECOND_Q += f'{SIGN} {coeff:.16f} {i+NUM}^ {i+NUM}' + '\n'

two_elec_dict = {}
for i in range(NUM):
    for j in range(NUM):
        for k in range(NUM):
            for l in range(NUM):
                coeff = TWO_ELEC_MO[i, j, k, l]/2
                if coeff < 1e-10:
                    continue
                SIGN = '+' if coeff > 0 else ''
                SECOND_Q += f'{SIGN} {coeff:.16f} {i}^ {k}^ {l} {j}' + '\n'
                SECOND_Q += f'{SIGN} {coeff:.16f} {i}^ {k+NUM}^ {l+NUM} {j}' + '\n'
                SECOND_Q += f'{SIGN} {coeff:.16f} {i+NUM}^ {k}^ {l} {j+NUM}' + '\n'
                SECOND_Q += f'{SIGN} {coeff:.16f} {i+NUM}^ {k+NUM}^ {l+NUM} {j+NUM}' + '\n'

pauli = JordanWignerMapper(SECOND_Q)

qc = QuantumCircuit(2*NUM)
for qubit in range(NUM-1):
    qc.x(qubit)
    qc.x(qubit+NUM)


for p_string, values in pauli.items():
    print(p_string, values)

    qc.barrier()

    for idx, p in p_string.items():
        if p.symbol == 'X':
            qc.h(idx)

        elif p.symbol == 'Y':
            qc.s(idx)
            qc.h(idx)

    qc.barrier()

    for i in range(2*NUM-1):
        qc.cx(i, i+1)

    qc.rz(values.real, 2*NUM-1)

    for i in reversed(range(2*NUM-1)):
        qc.cx(i, i+1)

    qc.barrier()

    for idx, p in p_string.items():
        if p.symbol == 'X':
            qc.h(idx)

        elif p.symbol == 'Y':
            qc.h(idx)
            qc.sdg(idx)




qc.draw('mpl')
