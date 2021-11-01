import math

import numpy as np

import vegas


def f(x):
    dx2 = sum(((_x - .5)**2 for _x in x))
    return math.exp(-dx2 * 100.) * 1013.2118364296088


integ = vegas.Integrator([[-1, 1], [0, 1], [0, 1], [0, 1]])

# step 1 -- adapt to f; discard results
integ(f, nitn=10, neval=2000)

# step 2 -- integ has adapted to f; keep results
result = integ(f, nitn=10, neval=10000)
print(result.summary())
print('result = %s    Q = %.2f' % (result, result.Q))
