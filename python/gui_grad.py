import math

import PySimpleGUI as sg
import matplotlib.pyplot as plt
from matplotlib import use as use_agg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

import numpy as np

import ctypes
from ctypes import *

from numpy.ctypeslib import ndpointer

# Use Tkinter Agg
use_agg("TkAgg")

# PySimplGUI window
layout = [
    [sg.Graph((640, 480), (0, 0), (640, 480), key="Graph")],
    [sg.Text("Calc type")],
    [
        sg.Text("Accurance"),
        sg.Input(key="accuracy", default_text="0.0001"),
        sg.Text("Iter count"),
        sg.Input(key="itercount", default_text="14"),
    ],
    [
        sg.Button("gradient descent method"),
        sg.Button("steepest descent method"),
        sg.Button("parabola"),
        sg.Button("goldenRatio"),
        sg.Button("brent"),
    ],
]

window = sg.Window("Matplotlib", layout, finalize=True)
fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, window["Graph"].Widget)
plot_widget = canvas.get_tk_widget()
plot_widget.grid(row=0, column=0)

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


def draw_gif(ax, function, label=""):
    if isinstance(function,list):
        dots = function
    else:
        dots = np.ctypeslib.as_array(
            (ctypes.c_double * 100).from_address(ctypes.addressof(function.contents))
        )
    n = 0
    print(dots)

    dots_y = []
    dots_x = []
    for i, d in enumerate(dots):
        if not i or dots[i] != dots[i - 1]:
            dots_y.append(func(d))
            dots_x.append(d)
            n += 1
        else:
            break

    ax.set_title(label)
    ax.plot(
        dots_x[:1], dots_y[:1], "*-", linewidth=0.6, markersize=4, label="opt", c="r"
    )
    ax.legend()

    def animate(i):
        ax.plot(
            dots_x[:i],
            dots_y[:i],
            "*-",
            linewidth=0.6,
            markersize=4,
            label="opt",
            c="r",
        )

    return animate, n


x = np.linspace(0, 5, 100).astype(float)
y = [func(i) for i in x]
i = n = 0

while True:

    event, values = window.read(timeout=200)
    if event == sg.WINDOW_CLOSED:
        break

    if event in funcs:
        accuracy = float(values["accuracy"])
        itercount = int(values["itercount"])
        f, n = draw_gif(ax, funcs[event](accuracy, itercount), event)
        i = 0

    # Reset ax
    ax.cla()
    ax.set_title("Sensor Data")
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.grid()
    ax.plot(x, y, c="b")

    if i < n:
        f(i)
        i += 1

    fig.canvas.draw()  # Draw curve really

window.close()
