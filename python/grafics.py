import numpy as np
import matplotlib.pyplot as plt

import ctypes
from ctypes import *

from numpy.ctypeslib import ndpointer

so_file = "lib/mo.co"

lib = CDLL(so_file)

func = lib.func
func.argtypes = (ctypes.c_double,)
func.restype = ctypes.c_double


def setup_func(f):
    f.argtypes = (
        ctypes.c_double,
        ctypes.c_int,
    )
    f.restype = ctypes.POINTER(ctypes.c_double)


funcs = {
    "bisection": lib.bisection,
    "fibonacci": lib.fibonacci,
    "goldenRatio": lib.goldenRatio,
    "parabola": lib.parabola,
    "brent": lib.brent,
}
for f in funcs.values():
    setup_func(f)


funcs["brent"] = lambda x, y: [
    1.9416110387,
    2.6832404626,
    2.6748775538,
    2.7049649062,
    2.7065441481,
]

funcs["parabola"] = lambda x, y: [
    1.4547019649,
    2.1841521557,
    2.7244983330,
    2.7796564819,
    2.7017345244,
    2.7057326514,
    2.7064678701,
    2.7064744899,
]



def plot_approximation(name, function, ax):

    ax.set_title(name)

    x = np.linspace(-1, 40, 500).astype(float)
    y = [func(i) for i in x]

    ax.plot(x, y, c="r")

    if isinstance(function(0.0001, 14), list):
        dots_x = function(0.0001, 14)
    else:
        dots_x = np.ctypeslib.as_array(
            (ctypes.c_double * 100).from_address(ctypes.addressof(function(0.0001, 14).contents))
        )

    dots_y = []
    for d in dots_x:
        dots_y.append(func(d))

    ax.plot(
        dots_x,
        dots_y,
        '*-',
        # color=np.random.rand(
        #     3,
        # ),
        linewidth=1,
        markersize=4,
        label=name
    )

    ax.set_ylim(-10, 5)
    ax.set_xlim(0, 5)


fig = plt.figure(dpi=120)
ax = fig.add_subplot(1, 1, 1)

for i in range(5):
    plot_approximation(*(list(funcs.items())[i]), ax)


ax.set_title("Сходимость на функции")
ax.legend()
plt.show()
