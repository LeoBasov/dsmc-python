from numba import cfunc, njit
import numba
import numpy.typing as npt

@cfunc("double(double[::3, ::2])")
def get_V(a):
    return (a[0][1] - a[0][0]) * (a[1][1] - a[1][0]) * (a[2][1] - a[2][0])

@njit((numba.int64[:], numba.int64, numba.int64), cache=True)
def swap(arr, pos1, pos2):
    temp = arr[pos2]
    arr[pos2] = arr[pos1]
    arr[pos1] = temp
    
@njit(cache=True)
def is_inside(position : npt.NDArray, box : npt.NDArray) -> bool:
    a : bool = position[0] >= box[0][0] and position[0] <= box[0][1]
    b : bool = position[1] >= box[1][0] and position[1] <= box[1][1]
    c : bool = position[2] >= box[2][0] and position[2] <= box[2][1]

    return a and b and c