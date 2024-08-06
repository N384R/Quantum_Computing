'''
Quantum Fourier Transform
'''

import numpy as np

def qft(qc, n):
    'Quantum Fourier Transform (QFT) circuit.'
    for i in range(n):
        qc.h(i)
        for j in range(i + 1, n):
            theta = 2 * np.pi / 2**(j - i + 1)
            qc.cp(theta, j, i)

def qft_dg(qc, n):
    'Inverse Quantum Fourier Transform (QFT^dagger) circuit'
    for i in reversed(range(n)):
        qc.h(i)
        for j in reversed(range(i)):
            theta = -2 * np.pi / 2**(i - j + 1)
            qc.cp(theta, i, j)
