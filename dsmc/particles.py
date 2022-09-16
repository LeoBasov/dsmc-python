import math
import numpy as np
from numba import njit

kb = 1.380649e-23

@njit
def box_muller(T):
    """
    Parameters
    ----------
    T : float
        temperature [K]

    Returns
    -------
    x : float
        m * v^2 / (2 * kb)
    """
    r1 = np.random.random()
    r2 = np.random.random()
    r3 = np.random.random()

    return T*(-np.log(r1) - np.log(r2) * math.pow(np.cos((np.pi/2.0)*r3), 2))

@njit
def x2velocity(x, mass):
    """
    Parameters
    ----------
    x : float
        m * v^2 / (2 * kb)
    mass : float
           particle mass [kg]

    Returns
    -------
    speed of particle : float
    """
    return math.sqrt((2.0/3.0) * x * kb /mass)

@njit
def get_vel(T, mass):
    """
    Parameters
    ----------
    T : float
        temperature [K]
    mass : float
           particle mass [kg]

    Returns
    -------
    velocity : np.array, shape = (3, 1)
    """
    return np.array([(-1)**(int(2*np.random.random())) * x2velocity(box_muller(T), mass) for _ in range(3)])

@njit
def get_velocities(T, mass, N):
    velocities = np.empty((N, 3), dtype=float)

    for i in range(N):
        velocities[i] = get_vel(T, mass)

    return velocities

class Particles:
    def __init__(self):
        self.velocities = None
        self.positions = None
        self.N = 0 # number of particles
