import numpy as np
import numpy.typing as npt
from numba import njit
from . import common as com

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
    a : bool = position[0] >= box[0][0] and position[0] <= box[0][1]
    b : bool = position[1] >= box[1][0] and position[1] <= box[1][1]
    c : bool = position[2] >= box[2][0] and position[2] <= box[2][1]

    return a and b and c
    
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
    runner = offset
    Nnew = 0
    for i in range(offset, offset + N):
        p = new_permutations[i]
        if _is_inside(positions[p], box):
            com.swap(new_permutations, i, runner)
            runner += 1
            Nnew += 1
            
    return new_permutations, Nnew

@njit
def _create_boxes(box):
    half = np.array([0.5*(box[i][0] + box[i][1]) for i in range(3)])
    
    child_geo1 = np.array(((half[0], box[0][1]), (half[1], box[1][1]), (half[2], box[2][1])))
    child_geo2 = np.array(((box[0][0], half[0]), (half[1], box[1][1]), (half[2], box[2][1])))
    child_geo3 = np.array(((box[0][0], half[0]), (box[1][0], half[1]), (half[2], box[2][1])))
    child_geo4 = np.array(((half[0], box[0][1]), (box[1][0], half[1]), (half[2], box[2][1])))
    
    child_geo5 = np.array(((half[0], box[0][1]), (half[1], box[1][1]), (box[2][0], half[2])))
    child_geo6 = np.array(((box[0][0], half[0]), (half[1], box[1][1]), (box[2][0], half[2])))
    child_geo7 = np.array(((box[0][0], half[0]), (box[1][0], half[1]), (box[2][0], half[2])))
    child_geo8 = np.array(((half[0], box[0][1]), (box[1][0], half[1]), (box[2][0], half[2])))
    
    return [child_geo1, child_geo2, child_geo3, child_geo4, child_geo5, child_geo6, child_geo7, child_geo8]

@njit
def _get_min_aspect_ratio(box, axis):
    half = np.array([0.5*(box[i][1] - box[i][0]) for i in range(3)])
    
    match axis:
        case 0:
            return min(half[0] / half[1], half[0] / half[2]);
        case 1:
            return min(half[1] / half[0], half[1] / half[2]);
        case 2:
            return min(half[2] / half[1], half[2] / half[0]);

@njit
def _devide(box, axis):
    half = 0.5*(box[axis][0] + box[axis][1])
    box1 = np.copy(box)
    box2 = np.copy(box)
    
    box1[axis][0] = box[axis][0]
    box1[axis][1] = half
    
    box2[axis][0] = half
    box2[axis][1] = box[axis][1]
    
    return (box1, box2)

@njit
def _create_combined_boxes(box, min_aspect_ratio):
    boxes = np.empty((15, 3, 2))
    boxes[0] = box
    N = 0
    Nold = 0
    q = 1
    
    for i in range(3):
        if _get_min_aspect_ratio(box, i) > min_aspect_ratio:
            for b in range(Nold, Nold + 2**N):
                new_boxes = _devide(boxes[b], i)
                boxes[q] = new_boxes[0]
                boxes[q + 1] = new_boxes[1]
                q += 2
            Nold += 2**N
            N += 1
            
    N = 2**N
    new_boxes = np.empty((N, 3, 2))
    
    for b in range(N):
        new_boxes[b] = boxes[Nold + b]
            
    return new_boxes
    
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
        self.min_aspect_ratio = 2.0/3.0
        
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
        if _is_resolved(self.cell_boxes[leaf_id], self.leafs[leaf_id].number_elements, self.w, self.sigma_T, self.Nmin, self.Nmax):
            
            self.leafs[leaf_id].id_first_child = self.cell_offsets[-1]
            
            new_boxes = _create_combined_boxes(self.cell_boxes[leaf_id], self.min_aspect_ratio)
            self.cell_offsets[-1] += len(new_boxes)
            self.leafs[leaf_id].number_children = len(new_boxes)
            
            for box in new_boxes:
                self.cell_boxes.append(box)
            
            #raise Exception()
        else:
            pass
           
        offset = 0 
           
        for i in range(self.leafs[leaf_id].number_children):
            new_leaf = Leaf()
            new_leaf.level = self.leafs[leaf_id].level + 1
            new_leaf.id_parent = leaf_id

            new_leaf.elem_offset = self.leafs[leaf_id].elem_offset + offset
            
            self.permutations, N = _sort(self.permutations, self.cell_boxes[self.leafs[leaf_id].id_first_child + i], positions, new_leaf.elem_offset, self.leafs[leaf_id].number_elements - offset)
            
            new_leaf.number_elements = N
            offset += N

            self.leafs.append(new_leaf)
