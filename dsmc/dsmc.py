import numpy as np
from numba import njit
from numba import prange
from . import particles as prt
from . import octree as oc
from . import common as com

@njit
def _push(velocities, positions, dt):
    old_positions = np.copy(positions)
    for p in prange(len(positions)):
        positions[p] = positions[p] + velocities[p]*dt
    return (velocities, positions, old_positions)

@njit
def _boundary(velocities, positions, old_positions, domain, boundary_conds):
    kept_parts = np.ones(positions.shape[0], dtype=np.uint)
    
    for p in prange(len(positions)):
        while not oc._is_inside(positions[p], domain) and kept_parts[p]:
            for i in range(3):
                if positions[p][i] < domain[i][0]:
                    if boundary_conds[i][0] == 0:
                        old_positions[p][i] = positions[p][i]
                        positions[p][i] = 2.0 * domain[i][0] - positions[p][i]
                        velocities[p][i] *= -1.0
                    elif boundary_conds[i][0] == 1 or boundary_conds[i][0] == 2:
                        kept_parts[p] = 0
                if positions[p][i] > domain[i][1]:
                    if boundary_conds[i][1] == 0:
                        old_positions[p][i] = positions[p][i]
                        positions[p][i] = 2.0 * domain[i][1] - positions[p][i]
                        velocities[p][i] *= -1.0
                    elif boundary_conds[i][1] == 1 or boundary_conds[i][0] == 2:
                        kept_parts[p] = 0
                        
    N = sum(kept_parts)
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
def _check_positions(velocities, positions, old_positions, domain):
    kept_parts = np.ones(positions.shape[0], dtype=np.uint)
    
    for i in prange(positions.shape[0]):
        if (not oc._is_inside(positions[i], domain)) and (not oc._is_inside(old_positions[i], domain)):
            kept_parts[i] = 0
    
    N = sum(kept_parts)
    p = 0
    new_velocities = np.empty((N, 3))
    new_positions = np.empty((N, 3))
    new_old_positions = np.empty((N, 3))
    
    for i in prange(positions.shape[0]):
        if kept_parts[i] == 1:
            new_velocities[p] = velocities[i]
            new_positions[p] = positions[i]
            new_old_positions[p] = old_positions[i]
            p += 1
        else:
            continue

    return (new_velocities, new_positions, new_old_positions)

@njit
def _check_created_particles(velocities, positions, obj):
    kept_parts = np.ones(positions.shape[0], dtype=np.uint)
    
    for i in prange(positions.shape[0]):
        if oc._is_inside(positions[i], obj):
            kept_parts[i] = 0
    
    N = sum(kept_parts)
    p = 0
    new_velocities = np.empty((N, 3))
    new_positions = np.empty((N, 3))
    
    for i in prange(positions.shape[0]):
        if kept_parts[i] == 1:
            new_velocities[p] = velocities[i]
            new_positions[p] = positions[i]
            p += 1
        else:
            continue

    return (new_velocities, new_positions)

@njit
def _object(velocities, positions, old_positions, coll_obj):
    for p in range(positions.shape[0]):
        if oc._is_inside(positions[p], coll_obj):
            for i in range(3):
                if (old_positions[p][i] < coll_obj[i][0]):
                    old_positions[p][i] = positions[p][i]
                    positions[p][i] = 2.0 * coll_obj[i][0] - positions[p][i]
                    velocities[p][i] *= -1.0
                
                if (old_positions[p][i] > coll_obj[i][1]):
                    old_positions[p][i] = positions[p][i]
                    positions[p][i] = 2.0 * coll_obj[i][1] - positions[p][i]
                    velocities[p][i] *= -1.0
                    
    return (velocities, positions, old_positions)

@njit
def _calc_prob(rel_vel : float, sigma_T : float, Vc : float, dt : float, w : float, N : int) -> np.single:
    """
    Parameters
    ----------
    vel1 : velocity
    vel2 : velocity
    sigma_T : float
        total cross section
    Vc : float
        cell volume
    w : float
        weight
    N : int
        number of particles

    Returns
    -------
    collision proability : float
    """
    return rel_vel * sigma_T * dt * w * N / Vc;

@njit
def _calc_post_col_vels(velocity1 : np.ndarray, velocity2 : np.ndarray, mass1 : float, mass2 : float, rel_vel_module : float, rand_number1 : float, rand_number2 : float) -> tuple[np.ndarray, np.ndarray]:
    mass12 = (mass1 + mass2)
    mass1_12 = (mass1 / mass12)
    mass2_12 = (mass2 / mass12)

    cos_xi = (2.0 * rand_number1 - 1.0)
    sin_xi = (np.sqrt(1.0 - cos_xi * cos_xi))
    epsilon = (2.0 * np.pi * rand_number2)

    centre_of_mass_velocity = (velocity1 * mass1 + velocity2 * mass2) * (1.0 / mass12)

    rel_velocity_new = np.empty((3, ))

    rel_velocity_new[0] = rel_vel_module * cos_xi
    rel_velocity_new[1] = rel_vel_module * sin_xi * np.cos(epsilon)
    rel_velocity_new[2] = rel_vel_module * sin_xi * np.sin(epsilon)

    return (centre_of_mass_velocity + rel_velocity_new * mass2_12 , centre_of_mass_velocity - rel_velocity_new * mass1_12)

@njit
def _update_velocities(permutations : np.ndarray, velocities : np.ndarray, mass : float, sigma_T : float, Vc : float, dt : float, w : float, offset : int, N : int) -> np.ndarray:
    for i in range(1, N, 2):
        p1 = permutations[offset + i - 1]
        p2 = permutations[offset + i]
        rel_vel = np.linalg.norm(velocities[p1] - velocities[p2])
        P = _calc_prob(rel_vel, sigma_T, Vc, dt, w, N)
        R = np.random.random(3)

        if R[0] < P:
            new_vels = _calc_post_col_vels(velocities[p1], velocities[p2], mass, mass, rel_vel, R[1], R[2])
            velocities[p1] = new_vels[0]
            velocities[p2] = new_vels[1]

    return velocities

@njit
def _update_vels(permutations : np.ndarray, velocities : np.ndarray, mass : float, sigma_T : float, dt : float, w : float, elem_offsets : np.ndarray, number_elements : np.ndarray, number_children : np.ndarray, cell_boxes : np.ndarray, Nleafs : int) -> np.ndarray:
    for i in range(Nleafs):
        if not number_children[i] and number_elements[i]:
            Vc = com.get_V(cell_boxes[i])
            velocities = _update_velocities(permutations, velocities, mass, sigma_T, Vc, dt, w, elem_offsets[i], number_elements[i])

    return velocities

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

class DSMC:
    def __init__(self):
        self.clear()

    def clear(self):
        self.particles = prt.Particles()
        self.octree = oc.Octree()
        self.w = None
        self.domain = None
        self.boundary_conds = np.array([[0, 0], [0, 0], [0, 0]], dtype=np.uint) # 0 = ela, 1 = open, 2 = inflow 
        self.boundary = Boundary()
        self.sigma_T = 3.631681e-19
        self.mass = None
        self.objects = []

    def advance(self, dt, collisions=True, octree=True):
        if self.domain is None:
            raise Exception("simulation domain not defined")
        if self.particles.N == 0:
            raise Exception("no particles in domain")
        if self.w == None:
            raise Exception("particle weight not set")
            
        for i in range(3):
            for j in range(2):
                if self.boundary_conds[i][j] == 2:
                    self.particles.inflow(self.mass, self.boundary.T[i][j], self.boundary.u[i][j], self.boundary.n[i][j], self.w, dt, self.domain, i, j)

        if octree:
            self.octree.build(self.particles.Pos)
        if collisions and octree:
            self.particles.VelPos = (self._update_velocities(dt), self.particles.Pos)
        velocities, positions, old_positions = _push(self.particles.Vel, self.particles.Pos, dt)
        velocities, positions, old_positions = _check_positions(velocities, positions, old_positions, self.domain)
        velocities, positions, old_positions = _boundary(velocities, positions, old_positions, self.domain, self.boundary_conds)
        
        for obj in self.objects:
            velocities, positions, old_positions  = _object(velocities, positions, old_positions, obj)
             
        self.particles.VelPos = (velocities, positions)

    def _update_velocities(self, dt):
        Nleafs : int = len(self.octree.leafs)
        elem_offsets : np.ndarray = np.array([leaf.elem_offset for leaf in self.octree.leafs], dtype=int)
        number_elements : np.ndarray = np.array([leaf.number_elements for leaf in self.octree.leafs], dtype=int)
        number_children : np.ndarray = np.array([leaf.number_children for leaf in self.octree.leafs], dtype=int)
        cell_boxes : np.ndarray = np.array([box for box in self.octree.cell_boxes])

        return _update_vels(self.octree.permutations, self.particles.Vel, self.mass, self.sigma_T, dt, self.w, elem_offsets, number_elements, number_children, cell_boxes, Nleafs)

    def create_particles(self, box, T, n, u = np.zeros(3)):
        box = np.array(box)
        N = int(round(com.get_V(box) * n / self.w))
        print("creating {} particles".format(N))
        self.particles.create_particles(box, self.mass, T, N, u)
        
        for obj in self.objects:
            self.particles.VelPos = _check_created_particles(self.particles.Vel, self.particles.Pos, obj)

        print("now containing {} particles, {} total".format(N, self.particles.N))

    def set_domain(self, domain):
        self.domain = np.array(domain)

    def set_mass(self, mass):
        self.mass = mass

    def set_weight(self, w):
        self.octree.w = w
        self.w = w
        
    def set_bc_type(self, boundary, bc_type):
        bound = _get_boundary(boundary)
        bc = _get_bc_type(bc_type)
        
        self.boundary_conds[bound[0]][bound[1]] = bc
        
        print("boundary [" + boundary + "] set to [" + bc_type + "]")
        
    def set_bc_values(self, boundary, T, n, u):
        i, j = _get_boundary(boundary)
        
        self.boundary.T[i][j] = T
        self.boundary.n[i][j] = n
        self.boundary.u[i][j] = u
        
        print("boundary [" + boundary + "] set to values T : {}, n : {}, u : {}".format(T, n, u))

    def add_object(self, coll_object):
        self.objects.append(np.array(coll_object))