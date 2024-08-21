from .fermion import Fermion
from .fermionic_sort import FermionicSort
from .pauli_operator import Pauli
from .pauli_string import PauliString, PauliStrings

class _JordanWigner:
    def __init__(self, fermionstring, maximum):
        self.fermion_string = [Fermion(fermion) for fermion in fermionstring]
        self.maximum = maximum
        self.pauli_strings = self.pauli_convert(self.fermion_string)

    def __iter__(self):
        return iter(self.pauli_strings)

    def _pauli_convert(self, fermion):
        paulistring1 = PauliString()
        paulistring2 = PauliString()
        z_string = PauliString()
        i_string = PauliString()

        num = fermion.num
        for i in range(num):
            z_string[i] = Pauli('Z')

        if num < self.maximum:
            for i in range(num+1, self.maximum+1):
                i_string[i] = Pauli('I')

        for key, value in z_string.items():
            paulistring1[key] = value
            paulistring2[key] = value

        paulistring1[num] = Pauli('X')
        operator_type = '-iY' if fermion.type == 'creation' else 'iY'
        paulistring2[num] = Pauli(operator_type)

        for key, value in i_string.items():
            paulistring1[key] = value
            paulistring2[key] = value

        return PauliStrings(paulistring1, paulistring2)

    def pauli_convert(self, fermionstring):
        result = None
        for fermion in fermionstring:
            pauli = self._pauli_convert(fermion)
            if result is None:
                result = pauli
            else:
                result *= pauli
        return result

    def __repr__(self):
        return f'{self.pauli_strings}'


class JordanWigner():
    def __init__(self, fermionstrings):
        self.fermion_strings = fermionstrings
        self.maximum = self._max_num(fermionstrings)
        compute_pauli = self.call_jordan_wigner(self.fermion_strings)
        self.pauli_strings = self._remove_zeros(compute_pauli)

    def __getitem__(self, key):
        return self.pauli_strings[key]

    def __iter__(self):
        return iter(self.pauli_strings)

    def call_jordan_wigner(self, fermionstrings):
        result = {}
        for fermionstring in fermionstrings:
            if fermionstring in ('-', '+'):
                sign = fermionstring
                continue

            coeff = float(fermionstring[0])/(2**(len(fermionstring[1:])))
            coeff = -coeff if sign == '-' else coeff
            pauli = _JordanWigner(fermionstring[1:], self.maximum)
            for p in pauli:
                c, p = self._pauli_arrange(p)
                if p in result:
                    result[p] += c * coeff
                else:
                    result[p] = c * coeff
        return result

    def _max_num(self, fermionstrings):
        maximum = 0
        for fermionstring in fermionstrings:
            if fermionstring in ('-', '+'):
                continue
            num = max(Fermion(fermion).num for fermion in fermionstring[1:])
            maximum = num if num > maximum else maximum
        return maximum

    def _pauli_arrange(self, pauli):
        coeff = 1
        for key, val in pauli.items():
            symbol = val.symbol
            count_m, count_i = symbol.count('-'), symbol.count('i')
            if count_m == 1:
                coeff *= -1
                pauli[key] *= Pauli('-I')
            if count_i == 1:
                coeff *= 1j
                pauli[key] *= Pauli('-iI')
        pauli.symbol = pauli.get_symbol(pauli)
        return coeff, pauli

    def _remove_zeros(self, pauli):
        for key, val in pauli.copy().items():
            if abs(val.real) < 1e-15 and abs(val.imag) < 1e-15:
                del pauli[key]
        # sorted_pauli = dict(sorted(pauli.items(), key=lambda x: str(x)))
        return pauli

    def __repr__(self):
        line = ''
        for key, values in self.pauli_strings.items():
            symbol = ''.join(f'{val}' for val in key.values())
            if values.real > 0 or values.imag > 0:
                line += '+ '
            else:
                line += '- '
            if values.real != 0 and values.imag != 0:
                line += f'{values}'
            elif values.real == 0 and values.imag != 0:
                line += f'({abs(values.imag):.16f})i'
            else:
                line += f'({abs(values.real):.16f})'
            line += ' ' + symbol + ' '
        return line

class JordanWignerMapper():
    def __init__(self, operator):
        self.fermi_sort = FermionicSort(operator)
        print(self.fermi_sort)
        self.pauli_strings = JordanWigner(self.fermi_sort).pauli_strings

    def __getitem__(self, key):
        return self.pauli_strings[key]

    def items(self):
        return self.pauli_strings.items()

    def __iter__(self):
        return iter(self.pauli_strings)

    def keys(self):
        return self.pauli_strings.keys()

    def values(self):
        return self.pauli_strings.values()

    def __mul__(self, other):
        return self.pauli_strings * other.pauli_strings

    def __add__(self, other):
        return self.pauli_strings + other.pauli_strings

    def __repr__(self):
        line = ''
        for key, values in self.pauli_strings.items():
            symbol = ''.join(f'{val}' for val in key.values())
            if values.real > 0 or values.imag > 0:
                line += '+ '
            else:
                line += '- '
            if values.real != 0 and values.imag != 0:
                line += f'{values}'
            elif values.real == 0 and values.imag != 0:
                line += f'({abs(values.imag):.15f})i'
            else:
                line += f'({abs(values.real):.15f})'
            line += ' ' + symbol + ' '
        return line

if __name__ == '__main__':
    jw = JordanWigner(['-', ['4.00', '1^']])
    print(jw)
