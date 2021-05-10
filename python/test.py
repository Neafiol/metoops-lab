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

function = func(0.005, 1, 3)
dots = np.ctypeslib.as_array(
    (ctypes.c_double * 4000).from_address(ctypes.addressof(function.contents))
)
old_val = 0
dots_x = []
dots_y = []
for x, y in dots.reshape((-1, 2)):
    if x + y == old_val or x + y == 0:
        break
    dots_x.append(x)
    dots_y.append(y)
    old_val = x + y
dots_x = np.array(dots_x)
dots_y = np.array(dots_y)

fig, ax = plt.subplots()
scale_x = max(abs(dots_x)) * 0.1
scale_y = max(abs(dots_y)) * 0.1
x, y = np.mgrid[min(dots_x)-scale_x:max(dots_x)+scale_x:0.1, min(dots_y)-scale_y:max(dots_y)+scale_y:0.1]

ax.contour(x, y, z3(x, y), levels=10)

ax.plot(
    dots_x,
    dots_y,
    "*-",
    linewidth=0.6,
    markersize=4,
    label="opt",
    c="r",
)

# fig.set_figwidth(6)  # ширина и
# fig.set_figheight(6)  # высота "Figure"

plt.show()
