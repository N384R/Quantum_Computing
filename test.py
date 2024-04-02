def combine_terms_with_coeffs(operator_list):
    terms_dict = {}  # terms를 키로, 계수의 합을 값으로 저장할 딕셔너리
    result = []

    for item in operator_list:
        if isinstance(item, list):  # 리스트 항목 처리
            coeff, *terms = item
            terms_tuple = tuple(terms)  # 리스트는 딕셔너리의 키가 될 수 없으므로 튜플로 변환

            # 동일한 terms에 대한 계수 합산
            if terms_tuple in terms_dict:
                terms_dict[terms_tuple] += coeff
            else:
                terms_dict[terms_tuple] = coeff
        else:  # 부호 항목 처리
            result.append(item)  # 부호는 결과 리스트에 직접 추가

    # 딕셔너리를 사용하여 최종 결과 리스트 구성
    for terms, coeff in terms_dict.items():
        if coeff != 0:  # 계수가 0이 아닌 항목만 결과에 추가
            result.append([coeff] + list(terms))

    return result

# 예제 입력
operator_list = ['+', [3, 'a', 'b'], '-', [2, 'c'], '+', [1, 'a', 'b'], '+', [4, 'd'], '-', [1, 'a', 'b']]

# 함수 호출 및 결과 출력
combined_terms = combine_terms_with_coeffs(operator_list)
print(combined_terms)