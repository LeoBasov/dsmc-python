from numba import cfunc

@cfunc("double(double[::3, ::2])")
def get_V(a):
    return (a[0][1] - a[0][0]) * (a[1][1] - a[1][0]) * (a[2][1] - a[2][0])