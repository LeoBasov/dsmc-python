import numpy as np
from numba import njit

@njit
def _check_if_parallel(v1, v2, diff=1e-6):
    V1 = np.copy(v1) / np.linalg.norm(v1)
    V2 = np.copy(v2) / np.linalg.norm(v2)
    
    return V1.dot(V2) > diff

@njit
def _intersect(l0, l1, p0, p1, p2):
    """
    Args:
        l0 : first position on line
        l1 : second position on line
        p0 : first position on plane
        p1 : first position on plane
        p2 : first position on plane
        
    Returns:
        (intersected, n_l, n_r, t)
    """
    n_l = l1 - l0
    n_p = np.cross((p1 - p0), (p2 - p1))
    
    if _check_if_parallel(n_l, n_p):
        return (False, n_l, n_p, 0.0)
    else:
        return (True, n_l, n_p, - ((l0 - p0).dot(n_p) / n_p.dot(n_l)))

@njit    
def _calc_nr(n_l, n_p):
    return n_l - 2.0 * (n_p.dot(n_l) / n_p.dot(n_p))*n_p

@njit     
def _reflect(vel, pos, pos_old, p0, p1, p2):
    intersected, n_l, n_p, t = _intersect(pos_old, pos, p0, p1, p2)
    
    if intersected:
        pos_old = pos_old + n_l*t
        n_r = _calc_nr(n_l, n_p)
        pos = pos_old  + (1.0 - t)*n_r
        vel = pos_old  + (np.linalg.norm(vel) / np.linalg.norm(n_r))*n_r

    return (vel, pos, pos_old)

@njit
def _get_plane(domain, i, j):
    if i == 0:
        p0 = np.array([domain[i][j], domain[1][0], domain[2][0]])
        p1 = np.array([domain[i][j], domain[1][1], domain[2][0]])
        p2 = np.array([domain[i][j], domain[1][0], domain[2][1]])
    elif i == 1:
        p0 = np.array([domain[0][0], domain[i][j], domain[2][0]])
        p1 = np.array([domain[0][1], domain[i][j], domain[2][0]])
        p2 = np.array([domain[0][0], domain[i][j], domain[2][1]])
    elif i == 2:
        p0 = np.array([domain[0][0], domain[1][0], domain[i][j]])
        p1 = np.array([domain[0][1], domain[1][0], domain[i][j]])
        p2 = np.array([domain[0][0], domain[1][1], domain[i][j]])
        
    return (p0, p1, p2)

@njit
def _boundary(velocities, positions, old_positions, domain, boundary_conds):
    kept_parts = np.ones(positions.shape[0], dtype=np.uint)
    
    for p in range(len(positions)):
        for i in range(3):
            for j in range(2):
                p0, p1, p2 = _get_plane(domain, i, j)
                if boundary_conds[i][j] == 0:
                    velocities[p], positions[p], old_positions[p] = _reflect(velocities[p], positions[p], old_positions[p], p0, p1, p2)
                elif boundary_conds[i][j] == 1 or boundary_conds[i][j] == 2:
                    if _intersect(old_positions[p], positions[p], p0, p1, p2)[0]:
                        kept_parts[p] = 0
                        
    N = int(sum(kept_parts))
    p = 0
    new_velocities = np.empty((N, 3))
    new_positions = np.empty((N, 3))
    new_old_positions = np.empty((N, 3))
    
    for i in range(positions.shape[0]):
        if kept_parts[i] == 1:
            new_velocities[p] = velocities[i]
            new_positions[p] = positions[i]
            new_old_positions[p] = old_positions[p]
            p += 1
        else:
            continue

    return (new_velocities, new_positions, new_old_positions)

@njit
def _get_boundary(boundary):
    if boundary == "xmin":
        return (0, 0)
    elif boundary == "xmax":
        return (0, 1)
    elif boundary == "ymin":
        return (1, 0)
    elif boundary == "ymax":
        return (1, 1)
    elif boundary == "zmin":
        return (2, 0)
    elif boundary == "zmax":
        return (2, 1)
    
@njit
def _get_bc_type(bc_type):
    if bc_type == "ela":
        return 0
    elif bc_type == "open":
        return 1
    elif bc_type == "inflow":
        return 2
    
class Boundary:
    def __init__(self):
        self.T = np.ones((3, 2))*300.0
        self.n = np.ones((3, 2))*1e+18
        self.u = np.zeros((3, 2, 3))
        self.boundary_conds = np.array([[0, 0], [0, 0], [0, 0]], dtype=np.uint) # 0 = ela, 1 = open, 2 = inflow 
        self.domain = None
        
    def boundary(self, velocities, positions, old_positions):
        return self._boundary(velocities, positions, old_positions, self.domain, self.boundary_conds)