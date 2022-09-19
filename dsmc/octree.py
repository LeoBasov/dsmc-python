import math
import numpy as np
from numba import njit

fmin = np.finfo(float).min
fmax = np.finfo(float).max

@njit
def _find_bounding_box(positions):
    box = np.array([[fmax, fmin], [fmax, fmin], [fmax, fmin]])
    
    for pos in positions:
        for i in range(3):
            if pos[i] < box[i][0]:
                box[i][0] = pos[i]
            if pos[i] > box[i][1]:
                box[i][1] = pos[i]
    
    return box

@njit    
def _calc_N_res(w, sigma_T, n):
    """
    Parameters
    ----------
    w : float
        particle weight
    sigma_t : float
              total cross section [m^2]
    n : float
        number density [1/m^3]
    """

    return int(round(np.sqrt(2.0) / (32.0 * w * sigma_T**3 * n**2)))

@njit  
def _calc_n(box, N, w):
    """
    Parameters
    ----------
    box : np.array(3, 3)
          cell
    N : int
        number of particles in cell
    w : float
        particle weight
        
    Returns
    -------
    number density : float
    """
    return np.prod(np.array([box[i][1] - box[i][0] for i in range(3)])) * N / w
    
class Leafs:
    def __init__(self):
        self.level = []
        self.elem_offset = []
        self.number_elements = []
        self.id_parent = []
        self.id_first_child = []
        self.number_children = []
        
    def add_leaf(self, level, elem_offset, number_elements, id_parent, id_first_child, number_children):
        self.level.append(level)
        self.elem_offset.append(elem_offset)
        self.number_elements.append(number_elements)
        self.id_parent.append(id_parent)
        self.id_first_child.append(id_first_child)
        self.number_children.append(number_children)
        
class Octree:
    def __init__(self):
        self.leafs = Leafs()
        self.permutations = None
        self.sigma_T = 3.631681e-19
        self.max_level = 20
        self.cell_bounding_boxes = []
        self.cell_offsets = []
        
    def build(self, positions):
        bounding_box = _find_bounding_box(positions)
        self.permutations = np.array([i for i in range(len(positions))])
        self._insert_first_leaf(positions, bounding_box)
        
        for i in range(self.max_level):
            pass
        
    def _insert_first_leaf(self, positions, box):
        level = 0;
        elem_offset = 0;
        number_elements = len(positions)
        id_parent = 0;
        id_first_child = 0;
        number_children = 0;

        self.cell_bounding_boxes.append(box)
        self.cell_offsets.append(0);
        self.cell_offsets.append(1);
        
        self.leafs.add_leaf(level, elem_offset, number_elements, id_parent, id_first_child, number_children)
        