import numpy as np
from . import particles as prt
from . import octree as oc
#from . import common as com

class Values:
    def __init__(self):
        self.number_elements = 0
        
    def calc_vals(self, positions, velocities, ids, box):
        #V = com.get_V(box)
        self.number_elements = len(ids)


def analyse_2d(positions, velocities, mesh2d, h):
    boxes = np.empty((len(mesh2d.cells), 3, 2))
    values = [Values() for _ in range(len(mesh2d.cells))]
    
    for y in range(mesh2d.n_cells2):
        for x in range(mesh2d.n_cells1):
            N = x + y * mesh2d.n_cells1
            
            boxes[N][0][0] = mesh2d.min1 + mesh2d.cell_size1*x
            boxes[N][0][1] = mesh2d.min1 + mesh2d.cell_size1*(x + 1)
            boxes[N][1][0] = mesh2d.min2 + mesh2d.cell_size2*y
            boxes[N][1][1] = mesh2d.min2 + mesh2d.cell_size2*(y + 1)
            boxes[N][2][0] = 0.0
            boxes[N][2][1] = h
            
            values[N].calc_vals(positions, velocities, mesh2d.cells[N], boxes[N])
        
    return (boxes, values)
                    

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
