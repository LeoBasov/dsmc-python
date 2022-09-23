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
    velocities = np.empty((N, 3))

    for i in range(N):
        velocities[i] = get_vel(T, mass)

    return velocities

@njit
def calc_temperature(velocities, mass):
    tot_e = 0.0

    if not len(velocities):
        return 0.0

    for i in range(len(velocities)):
        tot_e += 0.5 * mass * np.dot(velocities[i], velocities[i])

    return tot_e / ((3.0/2.0) * len(velocities) * kb)

@njit
def calc_positions(x, y, z, N):
    """
    Parameters
    ----------
    x : tuple(2), dtype = float
        xmin, xmax
    y : tuple(2), dtype = float
        ymin, ymax
    z : tuple(2), dtype = float
        zmin, zmax
    """
    positions = np.empty((N, 3))

    for i in range(N):
        positions[i][0] = x[0] + np.random.random() * (x[1] - x[0])
        positions[i][1] = y[0] + np.random.random() * (y[1] - y[0])
        positions[i][2] = z[0] + np.random.random() * (z[1] - z[0])

    return positions

class Particles:
    def __init__(self):
        self._velocities = None
        self._positions = None
        self._N = 0 # number of particles

    @property
    def Vel(self):
        return self._velocities

    @property
    def Pos(self):
        return self._positions

    @property
    def N(self):
        return self._N

    @property
    def VelPos(self):
        return (self._velocities, self._positions)

    @VelPos.setter
    def VelPos(self, vel_pos):
        assert len(vel_pos[0]) == len(vel_pos[1])

        self._velocities = vel_pos[0]
        self._positions = vel_pos[1]
        self._N = len(self._positions)
        
    def create_particles(self, X, mass, T, N):    
        if self._N == 0:
            self._velocities = get_velocities(T, mass, N)
            self._positions = calc_positions(X[0], X[1], X[2], N)
            self._N  = N
        else:
            self._velocities = np.concatenate((self._velocities, get_velocities(T, mass, N)))
            self._positions = np.concatenate((self._positions, calc_positions(X[0], X[1], X[2], N)))
            self._N  += N
