import ctypes
from ctypes import *

import matplotlib.pyplot as plt
import numpy as np

z1 = lambda x, y: 64 * x * x + 126 * x * y + 64 * y ** 2 - 10 * x + 30 * y + 13
z2 = lambda x, y: x * x + y * y - 2 * x + 5 * y
z3 = lambda x, y: 25 * x * x - 6 * x * y + y * y - 20 * x + 300 * y + 1

so_file = "lib/mog.co"

lib = CDLL(so_file)

func = lib.gradientDescent
func.argtypes = (ctypes.c_double,)
func.restype = ctypes.POINTER(ctypes.c_double)

function = func(0.005, 3, 3)
dots = np.ctypeslib.as_array(
    (ctypes.c_double * 4000).from_address(ctypes.addressof(function.contents))
)


# fig.set_figwidth(6)  # ширина и
# fig.set_figheight(6)  # высота "Figure"

plt.show()
