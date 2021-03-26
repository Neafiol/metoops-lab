import matplotlib
import numpy as np
import matplotlib.pyplot as plt

import ctypes
from ctypes import *

from matplotlib import animation
from matplotlib.lines import Line2D
from numpy.ctypeslib import ndpointer

so_file = "lib/mo.co"

lib = CDLL(so_file)

func = lib.func
func.argtypes = (ctypes.c_double,)
func.restype = ctypes.c_double

x = np.linspace(-1, 20, 100).astype(float)
y = [func(i) for i in x]


fig, ax = plt.subplots()
ax.plot(x, y, c="b")


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


def draw_gif(r, label=""):
    dots_x = np.ctypeslib.as_array(
        (ctypes.c_double * 100).from_address(ctypes.addressof(r.contents))
    )

    dots_y = []
    for d in dots_x:
        dots_y.append(func(d))

    ax.set_title(label)
    ax.plot(
        dots_x[:1], dots_y[:1], "*-", linewidth=0.6, markersize=4, label="opt", c="r"
    )
    ax.legend()

    def animate(i):
        # return ax.plot(x[:i], y[:i], c="r")
        (pl,) = ax.plot(
            dots_x[:i],
            dots_y[:i],
            "*-",
            linewidth=0.6,
            markersize=4,
            label="opt",
            c="r",
        )

        return (pl,)

    frames = []
    for i in range(len(dots_x)):
        for _ in range(10):
            frames.append(i)

    ani = animation.FuncAnimation(
        fig, animate, frames=frames, blit=True, interval=25, repeat=False
    )

    ani.save(f"line_{label}.gif", dpi=80, writer="imagemagick")


r = funcs["fibonacci"](0.0001, 14)

draw_gif(r, "fibonacci")
