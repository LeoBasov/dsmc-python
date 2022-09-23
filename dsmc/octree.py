import math
import numpy as np
import numpy.typing as npt
from numba import njit

fmin = np.finfo(float).min
fmax = np.finfo(float).max

@njit
def _find_bounding_box(positions : npt.NDArray) -> npt.NDArray:
    box = np.array([[fmax, fmin], [fmax, fmin], [fmax, fmin]])
    
    for pos in positions:
        for i in range(3):
            if pos[i] < box[i][0]:
                box[i][0] = pos[i]
            if pos[i] > box[i][1]:
                box[i][1] = pos[i]
    
    return box

@njit    
def _calc_N_res(w : float, sigma_T : float, n : float) -> int:
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
def _calc_n(box : npt.NDArray, N : float, w : float) -> float:
    """Calculates number density in cell
    
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

@njit  
def _is_resolved(box : npt.NDArray, N : int, w : float, sigma_T : float, Nmin : int, Nmax : int) -> bool:
    if N == 0:
        return False

    n = _calc_n(box, N, w)
    Nres = _calc_N_res(w, sigma_T, n)
    
    return N > 2 * min(Nmin, max(Nmin, Nres))

@njit
def _is_inside(position : npt.NDArray, box : npt.NDArray) -> bool:
    a : bool = position[0] > box[0][0] and position[0] <= box[0][1]
    b : bool = position[1] > box[1][0] and position[1] <= box[1][1]
    c : bool = position[2] > box[2][0] and position[2] <= box[2][1]

    return a and b and c

@njit    
def _swap(arr, pos1, pos2):
    temp = arr[pos2]
    arr[pos2] = arr[pos1]
    arr[pos1] = temp
    
@njit
def _sort(permutations : npt.NDArray, box : npt.NDArray, positions : npt.NDArray, offset : int, N : int) -> tuple[npt.NDArray, int]:
    '''sort particles in cell
    
    Parameters
    ----------
    box : np.ndarray
          cell
    positions : ndarray((3, ))
        particle positions
    offset : int
        number offset
    N : int
        number of particles to be considered
        
    Returns
    -------
    new_permutations : ndarray[int]
    N : int
        number of found positions
    '''
    new_permutations = np.copy(permutations)
    temp = np.empty((3,))
    runner = offset
    Nnew = 0
    for i in range(offset, offset + N):
        p = permutations[i]
        if _is_inside(positions[p], box):
            _swap(new_permutations, i, runner)
            runner += 1
            Nnew += 1
            
    return new_permutations, Nnew

@njit
def get_V(box):
    return (box[0][1] - box[0][0]) * (box[1][1] - box[1][0]) * (box[2][1] - box[2][0])
    
class Leaf:
    def __init__(self):
        self.level = 0
        self.elem_offset = 0
        self.number_elements = 0
        self.id_parent = None
        self.id_first_child = None
        self.number_children = 0
    
class Octree:
    def __init__(self):
        self.clear()
        
    def clear(self):
        self.cell_boxes = []
        self.leafs = []
        self.sigma_T = 3.631681e-19
        self.w = 1.0
        self.Nmin = 8
        self.Nmax = 64
        self.max_level = 10
        self.permutations = []
        self.cell_offsets = []
        self.level = 0
        
    def build(self, positions):
        self.clear()
        self._create_root(positions)
        self.permutations = np.array([i for i in range(len(positions))])
        
        for level in range(self.max_level):
            self.level += 1
            self.cell_offsets.append(self.cell_offsets[-1])
            for i in range(self.cell_offsets[level], self.cell_offsets[level + 1]):
                self._progress(i, positions)
   
            if self.cell_offsets[level + 1] == self.cell_offsets[level + 2]:
                break            
        
    def _create_root(self, positions):
        box = _find_bounding_box(positions)
        leaf = Leaf()
        leaf.number_elements = len(positions)
        
        self.cell_offsets += [0, 1]
        self.leafs.append(leaf)
        self.cell_boxes.append(box)
        
    def _progress(self, leaf_id, positions):
        leaf = self.leafs[leaf_id]
        if _is_resolved(self.cell_boxes[leaf_id], leaf.number_elements, self.w, self.sigma_T, self.Nmin, self.Nmax):
            leaf.number_children = 8
            leaf.id_first_child = leaf_id + 1
            self.cell_offsets[-1] += 8
            self._add_boxes(self.cell_boxes[leaf_id])
        else:
            pass
            
        offset = 0
            
        for i in range(leaf.number_children):
            new_leaf = Leaf()
            new_leaf.level = leaf.level + 1
            new_leaf.id_parent = leaf_id

            self.permutations, N = _sort(self.permutations, self.cell_boxes[leaf_id + 1 + i], positions, leaf.elem_offset, leaf.number_elements)
            
            new_leaf.number_elements = N
            new_leaf.elem_offset = leaf.elem_offset + offset
            offset += N

            self.leafs.append(new_leaf)
       
    def _add_boxes(self, box):
        half = np.array([0.5*(box[i][0] + box[i][1]) for i in range(3)])
        
        child_geo1 = np.array(((half[0], box[0][1]), (half[1], box[1][1]), (half[2], box[2][1])))
        child_geo2 = np.array(((box[0][0], half[0]), (half[1], box[1][1]), (half[2], box[2][1])))
        child_geo3 = np.array(((box[0][0], half[0]), (box[1][0], half[1]), (half[2], box[2][1])))
        child_geo4 = np.array(((half[0], box[0][1]), (box[1][0], half[1]), (half[2], box[2][1])))
        
        child_geo5 = np.array(((half[0], box[0][1]), (half[1], box[1][1]), (box[2][0], half[2])))
        child_geo6 = np.array(((box[0][0], half[0]), (half[1], box[1][1]), (box[2][0], half[2])))
        child_geo7 = np.array(((box[0][0], half[0]), (box[1][0], half[1]), (box[2][0], half[2])))
        child_geo8 = np.array(((half[0], box[0][1]), (box[1][0], half[1]), (box[2][0], half[2])))
        
        self.cell_boxes.append(child_geo1)
        self.cell_boxes.append(child_geo2)
        self.cell_boxes.append(child_geo3)
        self.cell_boxes.append(child_geo4)
        
        self.cell_boxes.append(child_geo5)
        self.cell_boxes.append(child_geo6)
        self.cell_boxes.append(child_geo7)
        self.cell_boxes.append(child_geo8)
