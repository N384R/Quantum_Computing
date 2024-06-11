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

    @staticmethod
    def generate_coeff(profile, coeff=1e-5):
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

    @staticmethod
    def mapping(profile, coeff):
        'Generate UCCSD ansatz Pauli strings'
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

    def ansatz(self, qc, profile, coeff):
        'Generate UCCSD ansatz circuit'
        chk = []
        for p_string, values in self.mapping(profile, coeff).items():
            chk.clear()
            for idx, p in p_string.items():
                if p.symbol == 'X':
                    qc.h(idx)
                elif p.symbol == 'Y':
                    qc.s(idx)
                    qc.h(idx)
                if p.symbol != 'I':
                    chk.append(idx)

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

class SpinFlipUCCSD(UCCSD):
    'Spin Flip Unitary Coupled Cluster Single and Double (UCCSD) ansatz'

    @staticmethod
    def generate_coeff(profile, coeff=1e-5):
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

    @staticmethod
    def mapping(profile, coeff):
        'Generate Spin Flip UCCSD ansatz Pauli strings'
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

class UCCGSD(UCCSD):
    'Unitary Coupled Cluster Generalized Single and Double (UCCGSD) ansatz'

    @staticmethod
    def generate_coeff(profile, coeff=1e-5):
        'Generate UCCGSD coefficients'
        no = profile.num_orb
        count = 0
        for i in range(no):
            for j in range(i, no):
                count += 2

        for i, j in product(range(no), range(no)):
            for _ in product(range(i, no), range(j, no)):
                count += 1

        for i in range(no-1):
            for _ in range(i+1, no):
                count += 1

        return [coeff] * count

    @staticmethod
    def mapping(profile, coeff):
        'Generate UCCGSD ansatz Pauli strings'
        uccsd_fermion = ''
        no = profile.num_orb
        value = iter(coeff)
        for i in range(no):
            for j in range(i+1, no):
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

        for i, j in product(range(no), range(no)):
            for k, l in product(range(i+1, no), range(j+1, no)):
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
