import math
import numpy
from numba import njit

kb = 1.380649e-23

@njit
def box_muller(T):
    r1 = np.random.random()
    r2 = np.random.random()
    r3 = np.random.random()

    return T*(-np.log(r1) - np.log(r2) * math.pow(np.cos((np.pi/2.0)*r3), 2))

@njit
def x2velocity(x, mass):
    return math.sqrt((2.0/3.0) * x * kb /mass)

@njit
def get_vel(T, mass):
    return np.array([(-1)**(int(2*np.random.random())) * x2velocity(box_muller(T), mass) for _ in range(3)])
