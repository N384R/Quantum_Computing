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

second_q = ''
for i in range(ORBITAL):
    for j in range(ORBITAL):
        second_q += f'{HCORE[i][j]} {i}^ {j} ' ## 더하기 부호 표시가 없음

for i in range(ORBITAL):
    for j in range(ORBITAL):
        for k in range(ORBITAL):
            for l in range(ORBITAL):
                second_q += f'{TWO_ELEC[i][j][k][l]} {i}^ {j}^ {k} {l} '

print(second_q)
pauli = JordanWignerMapper(second_q)
print(pauli)

#%%


qc = QuantumCircuit(2*ORBITAL)
for qubit in range(ORBITAL):
    qc.x(qubit)

qc.draw('mpl')
