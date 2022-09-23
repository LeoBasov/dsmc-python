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

class DSMC:
    def __init__(self):
        self.clear()
        
    def clear(self):
        self.particles = prt.Particles()
        self.octree = oc.Octree()
        self.w = None
        self.domain = None
        self.sigma_T = 3.631681e-19
        
    def advance(self, dt):
        if self.domain is None:
            raise Exception("simulation domain not defined")
        if self.particles.N == 0:
            raise Exception("no particles created")
        if self.w == None:
            raise Exception("particle weight not set")
            
        self.octree.build(self.particles.Pos)
        # update velocities
        positions = _push(self.particles.Vel, self.particles.Pos, dt)
        self.particles.VelPos = _boundary(self.particles.Vel, positions, self.domain)
        
    def create_particles(self, box, mass, T, n):
        N = int(round(n / self.w))
        print("creating {} particles".format(N))
        self.particles.create_particles(box, mass, T, N)
        
        print("now containing {} particles, {} total".format(N, self.particles.N))
        
    def set_domain(self, domain):
        self.domain = np.array(domain)
        
    def set_weight(self, w):
        self.octree.w = w
        self.w = w
