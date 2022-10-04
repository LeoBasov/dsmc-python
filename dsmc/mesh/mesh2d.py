from enum import Enum
import math
from numba import njit
import numpy as np

@njit
def _get_cell_id(val1, val2, n_cells1, n_cells2, min1, min2, cell_size):
    if (val1 < min1):
        return (False, 0)
    elif (val1 > (min1 + n_cells1 * cell_size)):
        return (False, 0)
    else:
        cell_id1 = math.floor((val1 - min1) / cell_size)

    if (val2 < min2):
        return (False, 0)
    elif (val2 > (min2 + n_cells2 * cell_size)):
        return (False, 0)
    else:
        cell_id2 = math.floor((val2 - min2) / cell_size)

    return (True, cell_id2 * n_cells1 + cell_id1)

@njit
def _sort(values1, values2, n_cells1, n_cells2, min1, min2, cell_size):
    inside = np.empty((len(values1)), np.bool_)
    ids = np.empty((len(values1)), np.int_)
    
    for i in range(len(values1)):
        inside[i], ids[i] = _get_cell_id(values1[i], values2[i], n_cells1, n_cells2, min1, min2, cell_size)
        
    return (inside, ids)


class Plane(Enum):
    XY = 0
    YZ = 1
    XZ = 2

class Mesh2d:
    def __init__(self):
        self.clear()
        
    def clear(self):
        self.cell_size = 1.0
        self.min1 = 0.0
        self.min2 = 0.0
        self.n_cells1 = 1
        self.n_cells2 = 1
        self.plane = Plane.XY
        
    def get_cell_id(self, position):
        match self.plane:
            case Plane.XY:
                return _get_cell_id(position[0], position[1], self.n_cells1, self.n_cells2, self.min1, self.min2, self.cell_size)
            case Plane.YZ:
                return _get_cell_id(position[1], position[2], self.n_cells1, self.n_cells2, self.min1, self.min2, self.cell_size)
            case Plane.XZ:
                return _get_cell_id(position[0], position[2], self.n_cells1, self.n_cells2, self.min1, self.min2, self.cell_size)