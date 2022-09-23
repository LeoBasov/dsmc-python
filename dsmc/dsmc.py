import math
import numpy as np
from numba import njit
from . import particles as prt
from . import octree as oc

@njit
def _push(velocities, positions, dt):  
    return positions + velocities*dt

@njit 
def _boundary(velocities, positions, domain):
    for p in range(len(positions)):
        while not oc._is_inside(positions[p], domain):
            for i in range(3):
                if positions[p][i] < domain[i][0]:
                    positions[p][i] = 2.0 * domain[i][0] - positions[p][i]
                    velocities[p][i] *= -1.0
                if positions[p][i] > domain[i][1]:
                    positions[p][i] = 2.0 * domain[i][1] - positions[p][i]
                    velocities[p][i] *= -1.0
                    
    return (velocities, positions)
 
@njit   
def _calc_prob(vel1 : np.ndarray, vel2 : np.ndarray, sigma_T : float, Vc : float, dt : float, w : float, N : int) -> np.single:
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
    return np.linalg.norm(vel1 - vel2) * sigma_T * dt * w * N / Vc;
  
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
    i = 1
    while i < N:
        p1 = permutations[offset + i - 1]
        p2 = permutations[offset + i]
        P = _calc_prob(velocities[p1], velocities[p2], sigma_T, Vc, dt, w, N)
        R = np.random.random(3)
        
        if R[0] < P:
            new_vels = _calc_post_col_vels(velocities[p1], velocities[p2], mass, mass, np.linalg.norm(velocities[p1] - velocities[p2]), R[1], R[2])
            velocities[p1] = new_vels[0]
            velocities[p2] = new_vels[1]
        
        i += 2
    
    return velocities
 
@njit   
def _update_vels(permutations : np.ndarray, velocities : np.ndarray, mass : float, sigma_T : float, dt : float, w : float, elem_offsets : np.ndarray, number_elements : np.ndarray, number_children : np.ndarray, cell_boxes : np.ndarray, Nleafs : int) -> np.ndarray:
    for i in range(Nleafs):
        if not number_children[i]:
            Vc = oc.get_V(cell_boxes[i])
            velocities = _update_velocities(permutations, velocities, mass, sigma_T, Vc, dt, w, elem_offsets[i], number_elements[i])
            
    return velocities

class DSMC:
    def __init__(self):
        self.clear()
        
    def clear(self):
        self.particles = prt.Particles()
        self.octree = oc.Octree()
        self.w = None
        self.domain = None
        self.sigma_T = 3.631681e-19
        self.mass = None
        
    def advance(self, dt):
        if self.domain is None:
            raise Exception("simulation domain not defined")
        if self.particles.N == 0:
            raise Exception("no particles created")
        if self.w == None:
            raise Exception("particle weight not set")
            
        self.octree.build(self.particles.Pos)
        
        self.particles.VelPos = (self._update_velocities(dt), self.particles.Pos)
        positions = _push(self.particles.Vel, self.particles.Pos, dt)
        self.particles.VelPos = _boundary(self.particles.Vel, positions, self.domain)
        
    def _update_velocities(self, dt):
        Nleafs : int = len(self.octree.leafs)
        elem_offsets : np.ndarray = np.array([leaf.elem_offset for leaf in self.octree.leafs], dtype=int)
        number_elements : np.ndarray = np.array([leaf.number_elements for leaf in self.octree.leafs], dtype=int)
        number_children : np.ndarray = np.array([leaf.number_children for leaf in self.octree.leafs], dtype=int)
        cell_boxes : np.ndarray = np.array([box for box in self.octree.cell_boxes])
        
        return _update_vels(self.octree.permutations, self.particles.Vel, self.mass, self.sigma_T, dt, self.w, elem_offsets, number_elements, number_children, cell_boxes, Nleafs)
        
    def create_particles(self, box, T, n, u = None):
        N = int(round(oc.get_V(box) * n / self.w))
        print("creating {} particles".format(N))
        self.particles.create_particles(np.array(box), self.mass, T, N)
        
        if u is not None:
            velocities = self.particles.Vel
            u = np.array(box)
            
            for i in range(len(velocities)):
                velocities[i] += u
            
            self.particles.VelPos = (velocities, self.particles.Pos)
        
        print("now containing {} particles, {} total".format(N, self.particles.N))
        
    def set_domain(self, domain):
        self.domain = np.array(domain)
        
    def set_mass(self, mass):
        self.mass = mass
        
    def set_weight(self, w):
        self.octree.w = w
        self.w = w
