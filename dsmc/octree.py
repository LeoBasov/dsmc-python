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
    n = _calc_n(box, N, w)
    N = _calc_N_res(w, sigma_T, n)
    
    return N > 2 * min(Nmin, max(Nmin, N))

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
def _sort(permutations : npt.NDArray, box : npt.NDArray, positions : npt.NDArray, offset : int, N : int) -> int:
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
