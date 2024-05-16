from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers import JordanWignerMapper

# 페르미온 연산자를 직접 정의합니다.
# 예를 들어, 간단한 2개의 페르미온 모드에 대한 헤밀토니안을 생각해봅시다.
# H = 1.0 * (a†0 a0) + 0.5 * (a†1 a1) + 0.3 * (a†0 a1) + 0.3 * (a†1 a0)
fermionic_op = FermionicOp({
    "+_0 -_0": -1.0,
    "+_1 -_1": -0.5,
    "+_2 -_2": -0.3,
    "+_0 +_0 -_0 -_0": 0.3,
    "+_1 +_1 -_1 -_1": 0.3,
    "+_0 +_1 -_0 -_1": 0.2,
    "+_1 +_2 -_2 -_1": 0.2
})

# JordanWignerMapper를 사용하여 페르미온 연산자를 파울리 연산자로 매핑합니다.
mapper = JordanWignerMapper()
qubit_op = mapper.map(fermionic_op)

# 결과 출력
print("Fermionic Operator:")
print(fermionic_op)
print("\nQubit Operator after Jordan-Wigner Mapping:")
print(qubit_op)
