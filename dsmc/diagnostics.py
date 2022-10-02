from numba import njit
import numpy as np
from . import particles as prt
from . import octree as oc

@njit
def sort_flat(positions, box, Nx, Ny):
    vals = np.zeros((Nx, Ny))
    dx = (box[0][1] - box[0][0]) / Nx
    dy = (box[1][1] - box[1][0]) / Ny
    
    for p in range(positions.shape[0]):
        pos = positions[p]
        X = box[0][0]
        Y = box[1][0]
        
        for x in range(Nx):
            for y in range(Ny):
                a = (pos[0] < (X + dx)) and (pos[0] >= X)
                b = (pos[1] < (Y + dy)) and (pos[1] >= Y)
                
                if a and b:
                    
                    vals[x][y] += 1
                
                Y += dy
            X += dx
    
    return vals.astype(np.int_)

def sort_bin(positions, axis, Nbin):
    bins = [[] for _ in range(Nbin)]
    b = 0
    box = oc._find_bounding_box(positions)
    dx = (box[axis][1] - box[axis][0]) / (Nbin - 1)
    xx = [dx]
    x = dx
    sub_pos = np.array([pos[axis] for pos in positions])
    sorted_pos = np.argsort(sub_pos)
    
    for i in range(len(sorted_pos)):
        p = sorted_pos[i]
        while positions[p][axis] > x:
            x += dx
            b += 1
            xx.append(x)
            
        bins[b].append(p)
        
    return bins, box, xx
    
def calc_n(bins, box, axis, w):
    Nbins = len(bins)
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
    
def calc_p(n, T):
    Nbins = len(n)
    p = np.empty((Nbins, ))
    
    for i in range(Nbins):
        p[i] = n[i]*T[i]*prt.kb
        
    return p
