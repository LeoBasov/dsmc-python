from numba import njit

@njit(float([[:] [:], [:]]))
def get_V(box):
    return (box[0][1] - box[0][0]) * (box[1][1] - box[1][0]) * (box[2][1] - box[2][0])