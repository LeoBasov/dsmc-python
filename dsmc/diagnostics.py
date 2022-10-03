import numpy as np
from numba import njit
from . import particles as prt
from . import octree as oc
from . import common as com

@njit
def _in_2d_box(pos, box):
    a = pos[0] >= box[0][0] and pos[0] <= box[0][1]
    b = pos[1] >= box[1][0] and pos[1] <= box[1][1]
    return a and b

@njit
def _set_up_2d_boxes(x0, x1, y0, y1, Nx, Ny):
    dx = (x1 - x0) / Nx
    dy = (y1 - y0) / Ny
    boxes = np.empty((Nx*Ny, 2, 2))
    y = y0
    
    for i in range(Ny):
        x = x0
        for j in range(Nx):
            k = j + i*Nx
            boxes[k][0][0] = x
            boxes[k][0][1] = x + dx
            boxes[k][1][0] = y
            boxes[k][1][1] = y + dy
            x += dx
        y += dy
    
    return boxes

@njit
def _sort_2d(positions, x0, x1, y0, y1, Nx, Ny):
    boxes = _set_up_2d_boxes(x0, x1, y0, y1, Nx, Ny) # [] : bix id, [][] : axis, [][][] : 0 min , 1 max
    permutations = [i for i in range(len(positions))]
    leafs = np.zeros((Nx * Ny, 2), dtype=np.int_) # leaf[0] : ofseat, leaf[1] : n-parts
    
    for i in range(len(leafs)):
        leafs[i][0] = leafs[i - 1][0] + leafs[i - 1][1] if i > 0 else 0
        leafs[i][1] = 0
        runner = leafs[i][0]
        for j in range(leafs[i][0], len(positions)):
            p = permutations[j]
            if _in_2d_box(positions[p], boxes[i]):
                com.swap(permutations, j, runner)
                runner += 1
                leafs[i][1] += 1

    return (permutations, leafs, boxes)
                    

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
