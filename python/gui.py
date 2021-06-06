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
        sg.Button("fibonacci"),
        sg.Button("bisection"),
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
    "parabola2": lib.parabola2,
    "brent2": lib.brent2,
}
for f in funcs.values():
    setup_func(f)


function = funcs["parabola2"](0.000001, 0)
dots = np.ctypeslib.as_array(
    (ctypes.c_double * 400).from_address(ctypes.addressof(function.contents))
)


def draw_gif(ax, function, label="", parabols=None):
    dots = np.ctypeslib.as_array(
        (ctypes.c_double * 100).from_address(ctypes.addressof(function.contents))
    )
    n = 0

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
        if parabols is not None:
            a, b, c = parabols[i]
            fparabol = lambda x: a * (x ** 2) + b * x + c

            px = np.linspace(0, 5, 100).astype(float)
            py = [fparabol(i) for i in x]

            ax.plot(px, py, c="y")

    return animate, n


x = np.linspace(0, 5, 100).astype(float)
y = [func(i) for i in x]
i = n = 0

while True:

    event, values = window.read(timeout=500)
    if event == sg.WINDOW_CLOSED:
        break

    if event in funcs:
        accuracy = float(values["accuracy"])
        itercount = int(values["itercount"])
        if event == "parabola":
            function = funcs["parabola2"](accuracy, itercount)
            dots = np.ctypeslib.as_array(
                (ctypes.c_double * 300).from_address(
                    ctypes.addressof(function.contents)
                )
            )
            parabols = dots.reshape((-1, 3))

        elif event == "brent":
            function = funcs["brent2"](accuracy, itercount)
            dots = np.ctypeslib.as_array(
                (ctypes.c_double * 300).from_address(
                    ctypes.addressof(function.contents)
                )
            )
            parabols = dots.reshape((-1, 3))
        else:
            parabols = None

        f, n = draw_gif(ax, funcs[event](accuracy, itercount), event, parabols)
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
