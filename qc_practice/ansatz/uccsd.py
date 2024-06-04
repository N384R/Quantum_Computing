'''
Unitary Coupled Cluster Single and Double (UCCSD) ansatz group

The class has the following methods:

generate_coeff: Generates UCCSD coefficients.
ansatz: Generates UCCSD ansatz circuit.
'''

from qc_practice.mapper.jordan_wigner import JordanWignerMapper

class UCCSD:
    'Unitary Coupled Cluster Single and Double (UCCSD) ansatz'
    def generate_coeff(self, profile, coeff=1e-5):
        'Generate UCCSD coefficients'
        n = profile.num_orb // 2
        return [coeff] * ((2 * n**2) + 2 * (n * (n - 1) // 2)**2 + n**4)

    def ansatz(self, profile, coeff):
        'Generate UCCSD ansatz circuit'
        uccsd_fermion = ''
        n = profile.num_orb
        value = iter(coeff)
        for i in range(n//2):
            for j in range(n//2, n):
                val = next(value)
                sign = '+' if val > 0 else '-'
                uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j} ' + '\n'
                sign = '-' if val > 0 else '+'
                uccsd_fermion += f'{sign} {abs(val):f} {j}^ {i} ' + '\n'

                val = next(value)
                sign = '+' if val > 0 else '-'
                uccsd_fermion += f'{sign} {abs(val):f} {i+n}^ {j+n} ' + '\n'
                sign = '-' if val > 0 else '+'
                uccsd_fermion += f'{sign} {abs(val):f} {j+n}^ {i+n} ' + '\n'

        for i in range(n//2):
            for j in range(n//2):
                for k in range(n//2, n):
                    for l in range(n//2, n):
                        val = next(value)
                        sign = '+' if val > 0 else '-'
                        uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j+n}^ {k} {l+n} ' + '\n'
                        sign = '-' if val > 0 else '+'
                        uccsd_fermion += f'{sign} {abs(val):f} {k}^ {l+n}^ {i} {j+n} ' + '\n'

        return JordanWignerMapper(uccsd_fermion)

class SpinFlipUCCSD:
    'Spin Flip Unitary Coupled Cluster Single and Double (UCCSD) ansatz'

    def generate_coeff(self, profile, coeff=1e-5):
        'Generate Spin Flip UCCSD coefficients'
        n = profile.num_orb // 2
        return [coeff] * ((2 * n**2) + 2 * (n * (n - 1) // 2)**2 + n**4 + n * (2 * n - 1))

    def ansatz(self, profile, coeff):
        'Generate Spin Flip UCCSD ansatz circuit'
        uccsd_fermion = ''
        n = profile.num_orb
        value = iter(coeff)
        for i in range(n//2):
            for j in range(n//2, n):
                val = next(value)
                sign = '+' if val > 0 else '-'
                uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j} ' + '\n'
                sign = '-' if val > 0 else '+'
                uccsd_fermion += f'{sign} {abs(val):f} {j}^ {i} ' + '\n'

                val = next(value)
                sign = '+' if val > 0 else '-'
                uccsd_fermion += f'{sign} {abs(val):f} {i+n}^ {j+n} ' + '\n'
                sign = '-' if val > 0 else '+'
                uccsd_fermion += f'{sign} {abs(val):f} {j+n}^ {i+n} ' + '\n'

        for i in range(n//2):
            for j in range(n//2):
                for k in range(n//2, n):
                    for l in range(n//2, n):
                        val = next(value)
                        sign = '+' if val > 0 else '-'
                        uccsd_fermion += f'{sign} {abs(val):f} {i}^ {j+n}^ {k} {l+n} ' + '\n'
                        sign = '-' if val > 0 else '+'
                        uccsd_fermion += f'{sign} {abs(val):f} {k}^ {l+n}^ {i} {j+n} ' + '\n'

                        if i == j and k == l:
                            val = next(value)
                            sign = '+' if val > 0 else '-'
                            uccsd_fermion += f'{sign} {abs(val):f} {i}^ {k+n}^ {k} {i+n} ' + '\n'
                            sign = '-' if val > 0 else '+'
                            uccsd_fermion += f'{sign} {abs(val):f} {k}^ {i+n}^ {i} {k+n} ' + '\n'

        return JordanWignerMapper(uccsd_fermion)
