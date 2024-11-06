from concurrent.futures import ProcessPoolExecutor
import scipy.optimize as opt
import time

def summation(val):
    'Summation function'
    time.sleep(0.1)
    return val

def polynomial(coeff):
    'Polynomial function'
    x = sum(summation(val) for val in coeff)
    print(f'x value: {x}', end='\r', flush=True)
    return x**2 - 4*x + 4

def mp_polynomial(coeff):
    'Polynomial function'
    with ProcessPoolExecutor(4) as pool:
        result = list(pool.map(summation, coeff))
    x = sum(result)
    print(f'x value: {x}', end='\r', flush=True)
    return x**2 - 4*x + 4

def call_optimizer(func, coeff):
    'Optimize the coefficients'
    start_time = time.time()
    result = opt.minimize(func, coeff, method='COBYLA')

    print(f'minimum: {result.fun:.2f} x: {result.x}')
    print(f'time passed: {time.time() - start_time:.2f}')

def call_twice(func, coeff):
    'Call optimizer twice'
    call_optimizer(func, coeff)
    call_optimizer(func, coeff)

# print('\nSingle process')
# call_optimizer(polynomial, [1, 2, 3, 4, 5])

# print('\nMultiprocess')
# call_optimizer(mp_polynomial, [1, 2, 3, 4, 5])

print('\nCall twice')
call_twice(mp_polynomial, [1, 2, 3, 4, 5])
