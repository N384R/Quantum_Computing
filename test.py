import random

def generate_block():
    # 상수는 -999.999에서 999.999 범위 내에서 랜덤하게 생성
    constant = round(random.uniform(-999.999, 999.999), 3)
    
    # 연산자 수는 1에서 10 사이로 랜덤하게 결정
    operators_count = random.randint(1, 10)
    
    # 연산자 생성: 0부터 9 사이의 숫자를 연산자로 하며, 각 연산자는 소멸 연산자 또는 생성 연산자일 수 있음
    operators = [f"{random.randint(0, 9)}^" if random.random() < 0.5 else str(random.randint(0, 9)) for _ in range(operators_count)]
    
    # 블럭을 문자열로 조합
    block = f"{constant} " + " ".join(operators)
    
    return block

# 100개의 블럭을 생성하고 라인으로 조합
blocks = [generate_block() for _ in range(100)]

# 첫 번째 블록 다음부터는 부호로 시작하도록 조정
line = blocks[0] + " " + " - ".join(blocks[1:])

line[:500]  # 결과 확인을 위해 처음 500자만 출력

def replace_double_negatives(line):
    return line.replace("- -", "+ ")

# 치환 함수 적용
corrected_line = replace_double_negatives(line)

# 결과의 일부 출력
print(corrected_line)