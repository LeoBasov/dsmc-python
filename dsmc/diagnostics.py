import numpy as np
from numba import njit
from . import particles as prt
from . import octree as oc

def sort_bin(positions, axis, Nbin):
    bins = [[] for _ in range(Nbin)]
    b = 0
    box = oc._find_bounding_box(positions)
    dx = (box[axis][1] - box[axis][0]) / (Nbin - 1)
    x = dx
    sub_pos = np.array([pos[axis] for pos in positions])
    sorted_pos = np.argsort(sub_pos)
    
    for i in range(len(sorted_pos)):
        p = sorted_pos[i]
        while positions[p][axis] > x:
            x += dx
            b += 1
            
        bins[b].append(p)
        
    return bins, box
    
def calc_n(bins, box, axis, w):
    Nbins = len(bins)
    dx = (box[axis][1] - box[axis][0]) / Nbins
    V = 1
    n = np.empty((Nbins, ))
    for i in range(3):
        if i == axis:
            V *= (box[i][1] - box[i][0]) / Nbins
        else:
            V *= (box[i][1] - box[i][0])
            
    for i in range(Nbins):
        n[i] = len(bins[i]) * w / V
        
    return n
    
def calc_T(bins, velocities, mass):
    Nbins = len(bins)
    T = np.empty((Nbins, ))
    
    for i in range(Nbins):
        vels = np.array([velocities[p] for p in bins[i]])
        T[i] = prt.calc_temperature(vels, mass)
        
    return T
