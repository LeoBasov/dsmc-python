import numpy as np
from numba import njit

@njit
def _check_if_parallel(v1, v2, diff=1e-6):
    V1 = np.copy(v1) / np.linalg.norm(v1)
    V2 = np.copy(v2) / np.linalg.norm(v2)
    
    return V1.dot(V2) > diff

def _intersect(l0, l1, p0, p1, p2):
    """
    Args:
        l0 : first position on line
        l1 : second position on line
        p0 : first position on plane
        p1 : first position on plane
        p2 : first position on plane
    """
    n_l = l1 - l0
    n_p = (p1 - p0).cross(p2 - p1)