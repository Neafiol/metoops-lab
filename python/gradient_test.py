import ctypes

class EmxArray(ctypes.Structure):
    """ creates a struct to match emxArray_real_T """

    _fields_ = [('data', ctypes.POINTER(ctypes.c_double)),
                ('size', ctypes.POINTER(ctypes.c_int)),
                ('allocatedSize', ctypes.c_int),
                ('numDimensions', ctypes.c_int),
                ('canFreeData', ctypes.c_bool)]

data = (1.3, 3.5, 2.7, 4.1)
L = len(data)

e = EmxArray()
e.data = (ctypes.c_double * L)(*data)
e.size = (ctypes.c_int * 1)(L)