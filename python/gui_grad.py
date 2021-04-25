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

functions = {
    "64x^2 + 126xy + 64y^2 - 10x + 30y + 13": {
        "func": lambda x, y: 64 * x * x
        + 126 * x * y
        + 64 * y ** 2
        - 10 * x
        + 30 * y
        + 13,
        "c_func_num": 1,
    },
    "x^2 + y^2 - 2x + 5y": {
        "func": lambda x, y: x * x + y * y - 2 * x + 5 * y,
        "c_func_num": 2,
    },
    "25x^2 - 6xy + y^2 - 20x + 300y + 1": {
        "func": lambda x, y: 25 * x * x - 6 * x * y + y * y - 20 * x + 300 * y + 1,
        "c_func_num": 3,
    },
}


layout = [
    [sg.Graph((640, 480), (0, 0), (640, 480), key="Graph")],
    [
        sg.Listbox(
            list(functions.keys()),
            key="func",
            select_mode=sg.LISTBOX_SELECT_MODE_SINGLE,
            default_values=[list(functions.keys())[0]],
            size=(35, 3),
            bind_return_key=True,
            auto_size_text=True,
            enable_events=True,
        ),
        sg.Button("metod 1"),
        sg.Button("metod 2"),
        sg.Button("metod 3"),
        sg.Text("Accurance"),
        sg.Input(key="accuracy", default_text="0.0001", size=(10, 1)),
    ],
]

window = sg.Window("Matplotlib", layout, finalize=True)
fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, window["Graph"].Widget)
plot_widget = canvas.get_tk_widget()
plot_widget.grid(row=0, column=0)


lib = CDLL("lib/mog.co")
get_gradient_trek = lib.gradientDescent
get_gradient_trek.argtypes = (ctypes.c_double,)
get_gradient_trek.restype = ctypes.POINTER(ctypes.c_double)

function = get_gradient_trek(0.005, 3, 3)
dots = np.ctypeslib.as_array(
    (ctypes.c_double * 4000).from_address(ctypes.addressof(function.contents))
)


def draw_gif(ax, function_data, label="", accuracy=0.005, metod_num=1):
    function = get_gradient_trek(accuracy, metod_num, function_data["c_func_num"])
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

    scale_x = max(abs(dots_x)) * 0.1
    scale_y = max(abs(dots_y)) * 0.1

    x, y = np.mgrid[
        min(dots_x) - scale_x : max(dots_x) + scale_x : 0.1,
        min(dots_y) - scale_y : max(dots_y) + scale_y : 0.1,
    ]

    ax.set_title(label)
    ax.legend()

    def animate(i):
        ax.contour(x, y, function_data["func"](x, y), levels=10)
        ax.plot(
            dots_x[:i],
            dots_y[:i],
            "*-",
            linewidth=0.6,
            markersize=4,
            label="opt",
            c="r",
        )

    return animate, len(dots_x)


i = n = 0
while True:

    event, values = window.read(timeout=100)
    if event == sg.WINDOW_CLOSED:
        break

    if event in ["metod 1", "metod 2", "metod 3"]:
        accuracy = float(values["accuracy"])
        func_name = values["func"][0]
        if event == "metod 1":
            draw_func, n = draw_gif(ax, functions[func_name], func_name, accuracy, 1)
        elif event == "metod 2":
            draw_func, n = draw_gif(ax, functions[func_name], func_name, accuracy, 2)
        else:
            draw_func, n = draw_gif(ax, functions[func_name], func_name, accuracy, 3)

        i = 0

    ax.cla()
    ax.set_title("Sensor Data")
    ax.set_xlabel("X axis")
    ax.set_ylabel("Y axis")
    ax.grid()

    if i < n:
        draw_func(i)
        i += 1

    fig.canvas.draw()  # Draw curve really

window.close()
