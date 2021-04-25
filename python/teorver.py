import random

import numpy
import numpy as np


def e_plus_d(epa, pow):
    ex = 0
    dx = 0
    for e, p in zip(epa, pow):
        ex += e * p

    for e, p in zip(epa, pow):
        dx = p * (e - ex) ** 2

    print(ex, dx)
    print(ex + dx)


# e_plus_d([-4, -2, 0, 3, 5], [0.1, 0.3, 0.3, 0.2, 0.1])


def get_p(x):
    if x < 0.5:
        return 2 * x
    return -2 / 3 * x + 4 / 3


def get_p3(x):
    return x


def get_attemps_count():
    n = 0
    while n<5 and random.random() < 1/(5 - n):
        n += 1
    return n


epa = []
pow = []
N = 10000000
for i in range(N):
    # x = random.random()*6 - 2
    epa.append(get_attemps_count())
    # pow.append(1 / N)

results = np.array(epa)
# print(results)
a = results.mean()
b = np.var(results)
print(a, b)
print(a + b)
