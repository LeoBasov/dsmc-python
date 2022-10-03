from numba import cfunc, njit

@cfunc("double(double[::3, ::2])")
def get_V(a):
    return (a[0][1] - a[0][0]) * (a[1][1] - a[1][0]) * (a[2][1] - a[2][0])

@njit 
def swap(arr, pos1, pos2):
    temp = arr[pos2]
    arr[pos2] = arr[pos1]
    arr[pos1] = temp