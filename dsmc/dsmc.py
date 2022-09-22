import math
import numpy as np
from numba import njit
from . import particles as prt
from . import octree as oc

@njit
def _push(velocities, positions, dt):
    return positions + velocities*dt

class DSMC:
    def __init__(self):
        self.clear()
        
    def clear(self):
        self.particles = prt.Particles()
        self.octree = oc.Octree()
        self.w = None
        self.domain = None
        
    def advance(self, dt):
        if self.domain is None:
            raise Exception("simulation domain not defined")
        if self.particles.N == 0:
            raise Exception("no particles created")
        if self.w == None:
            raise Exception("particle weight not set")
            
        self.octree.build(self.particles.Pos)
        # update velocities
        self.particles.VelPos = (self.particles.Vel, _push(self.particles.Vel, self.particles.Pos, dt))
        # do boundary
        
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
