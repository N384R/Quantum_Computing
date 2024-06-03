'''
Unitary Coupled Cluster Single and Double (UCCSD) ansatz group
'''

from qc_practice.mapper.jordan_wigner import JordanWignerMapper

class UCCSD:
    '''
    Unitary Coupled Cluster Single and Double (UCCSD) ansatz
    '''
    def ansatz(self, coeff, profile):
        '''
        Generate UCCSD ansatz circuit
        '''
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
    '''
    Spin Flip Unitary Coupled Cluster Single and Double (UCCSD) ansatz
    '''
    def ansatz(self, coeff, profile):
        '''
        Generate Spin Flip UCCSD ansatz circuit
        '''
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
