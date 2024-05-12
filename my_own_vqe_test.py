#%%
from pyscf import gto, scf
from qc_practice import JordanWignerMapper
from qiskit import QuantumCircuit

mol = gto.M(atom = 'H 0 0 0; H 0 0 0.74', basis = 'sto-3g')
rhf = scf.RHF(mol)
EN = rhf.kernel()
HCORE = rhf.get_hcore()
TWO_ELEC = mol.intor('int2e')
OCC = rhf.get_occ()
ORBITAL = HCORE.shape[0]

SECOND_Q = ''
for i in range(ORBITAL):
    for j in range(ORBITAL):
        SIGN = ' +' if HCORE[i][j] > 0 else ''
        SECOND_Q += f'{SIGN} {HCORE[i][j]} {i}^ {j} '

for i in range(ORBITAL):
    for j in range(ORBITAL):
        for k in range(ORBITAL):
            for l in range(ORBITAL):
                SIGN = ' +' if TWO_ELEC[i][j][k][l] > 0 else ''
                SECOND_Q += f'{SIGN} {TWO_ELEC[i][j][k][l]} {i}^ {j}^ {k} {l} '

pauli = JordanWignerMapper(SECOND_Q)
print(pauli)

#%%


# qc = QuantumCircuit(2*ORBITAL)
# for qubit in range(ORBITAL):
#     qc.x(qubit)

# qc.draw('mpl')
