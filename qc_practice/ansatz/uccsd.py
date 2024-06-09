'''
Unitary Coupled Cluster Single and Double (UCCSD) ansatz group

The class has the following methods:

generate_coeff: Generates UCCSD coefficients.
ansatz: Generates UCCSD ansatz circuit.
'''
from itertools import product
from qc_practice.mapper.jordan_wigner import JordanWignerMapper

class UCCSD:
    'Unitary Coupled Cluster Single and Double (UCCSD) ansatz'
    def generate_coeff(self, profile, coeff=1e-5):
        'Generate UCCSD coefficients'
        no = profile.num_orb
        ne = profile.num_elec
        count = 0
        for _ in product(range(ne//2), range(ne//2, no)):
            count += 2
        for _ in product(range(ne//2), range(ne//2),
                        range(ne//2, no), range(ne//2, no)):
            count += 1

        return [coeff] * count

    def ansatz(self, profile, coeff):
        'Generate UCCSD ansatz circuit'
        uccsd_fermion = ''
        no = profile.num_orb
        ne = profile.num_elec
        value = iter(coeff)
        for i, j in product(range(ne//2), range(ne//2, no)):
            val = next(value)
            sign = '+' if val > 0 else '-'
            uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j} ' + '\n'
            sign = '-' if val > 0 else '+'
            uccsd_fermion += f'{sign} {abs(val):f} {j}^ {i} ' + '\n'

            val = next(value)
            sign = '+' if val > 0 else '-'
            uccsd_fermion += f'{sign} {abs(val):f} {i+no}^ {j+no} ' + '\n'
            sign = '-' if val > 0 else '+'
            uccsd_fermion += f'{sign} {abs(val):f} {j+no}^ {i+no} ' + '\n'

        for i, j, k, l in product(range(ne//2), range(ne//2),
                                  range(ne//2, no), range(ne//2, no)):
            val = next(value)
            sign = '+' if val > 0 else '-'
            uccsd_fermion += f'{sign} {abs(val):f} {k}^ {l+no}^ {i} {j+no} ' + '\n'
            sign = '-' if val > 0 else '+'
            uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j+no}^ {k} {l+no} ' + '\n'

        return JordanWignerMapper(uccsd_fermion)

class SpinFlipUCCSD:
    'Spin Flip Unitary Coupled Cluster Single and Double (UCCSD) ansatz'

    def generate_coeff(self, profile, coeff=1e-5):
        'Generate Spin Flip UCCSD coefficients'
        no = profile.num_orb
        ne = profile.num_elec
        count = 0
        for _ in product(range(ne//2), range(ne//2, no)):
            count += 2
        for _ in product(range(ne//2), range(ne//2),
                                  range(ne//2, no), range(ne//2, no)):
            count += 1
        for i in range(no-1):
            for _ in range(i+1, no):
                count += 1

        return [coeff] * count

    def ansatz(self, profile, coeff):
        'Generate Spin Flip UCCSD ansatz circuit'
        uccsd_fermion = ''
        no = profile.num_orb
        ne = profile.num_elec
        value = iter(coeff)
        for i in range(ne//2):
            for j in range(ne//2, no):
                val = next(value)
                sign = '+' if val > 0 else '-'
                uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j} ' + '\n'
                sign = '-' if val > 0 else '+'
                uccsd_fermion += f'{sign} {abs(val):f} {j}^ {i} ' + '\n'

                val = next(value)
                sign = '+' if val > 0 else '-'
                uccsd_fermion += f'{sign} {abs(val):f} {i+no}^ {j+no} ' + '\n'
                sign = '-' if val > 0 else '+'
                uccsd_fermion += f'{sign} {abs(val):f} {j+no}^ {i+no} ' + '\n'

        for i, j, k, l in product(range(ne//2), range(ne//2),
                                  range(ne//2, no), range(ne//2, no)):
            val = next(value)
            sign = '+' if val > 0 else '-'
            uccsd_fermion += f'{sign} {abs(val):f} {k}^ {l+no}^ {i} {j+no} ' + '\n'
            sign = '-' if val > 0 else '+'
            uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j+no}^ {k} {l+no} ' + '\n'

        for i in range(no-1):
            for j in range(i+1, no):
                val = next(value)
                sign = '+' if val > 0 else '-'
                uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j+no}^ {j} {i+no} ' + '\n'
                sign = '-' if val > 0 else '+'
                uccsd_fermion += f'{sign} {abs(val):f} {j}^ {i+no}^ {i} {j+no} ' + '\n'

        return JordanWignerMapper(uccsd_fermion)
