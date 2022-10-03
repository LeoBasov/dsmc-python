import numpy as np
from . import particles as prt
from . import octree as oc

def _sort_2d(positions, x0, x1, y0, y1, Nx, Ny):
    dx = (x1 - x0) / Nx
    dy = (y1 - y0) / Ny
    sorted_vals = [[[] for _ in range(Ny)] for _ in range(Nx)]
    
    for p in range(len(positions)):
        pos = positions[p]
        x = x0
        y = y0
        brk = False
        
        for i in range(Nx):
            if brk:
                break
            y = y0
            
            for j in range(Ny):
                a = pos[0] >= x and pos[0] <= (x + dx)
                b = pos[1] >= y and pos[1] <= (y + dy)
                if a and b:
                    sorted_vals[i][j].append(p)
                    brk = True
                    break
                y += dy
            x += dx

    return sorted_vals
                    

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
