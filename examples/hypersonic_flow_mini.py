import dsmc
import dsmc.writer as wrt
import dsmc.octree as oc
import time
import numpy as np

if __name__ == '__main__':
    # general parameters
    solver = dsmc.DSMC()
    domain = [(-3.0, 3.0), (-1.5, 1.5), (-0.025, 0.025)]
    obj = [(-0.25, 0.25), (-0.25, 0.25), (-0.5, 0.5)]
    dt = 1e-6
    w = 0.2 * 2.3e+15
    mass = 6.6422e-26
    T =  273.0
    n = 2.6e+19
    u = np.array([0.0, -3043.0, 0.0])
    niter = 500
    
    h = domain[2][1] - domain[2][0]
    
    # setup solver
    solver.set_domain(domain)
    solver.set_weight(w)
    solver.set_mass(mass)
    
    # set object
    solver.add_object(obj)
    
    # trees
    tree_inner = oc.Octree()
    tree_outer = oc.Octree()
    
    inner_pos = np.array([[-0.25, -0.25, -0.025], [0.25, 0.25, 0.025]])
    outer_pos = np.array([[-3.0, -1.5, -0.025], [3.0, 1.5, 0.025]])
    
    tree_inner.build(inner_pos)
    tree_outer.build(outer_pos)
    
    # create particles
    positions = np.zeros((2, 3))
    velocities = np.zeros((2, 3))
    
    positions[0][0] = -0.01
    positions[0][1] = 1.4
    positions[0][2] = -0.01
    
    positions[1][0] = 0.01
    positions[1][1] = 1.3
    positions[1][2] = 0.01
    
    velocities[0][1] = -3043.0
    velocities[1][1] = -3043.0
    
    solver.particles.VelPos = (velocities, positions)
    
    # start timing
    start_time = time.time()
    
    for it in range(niter):
        print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)
        solver.advance(dt)
        wrt.write_positions(solver.particles, "pos_{}.csv".format(it))
        
    wrt.write_buttom_leafs(tree_inner, "innter.vtu")
    wrt.write_buttom_leafs(tree_outer, "outer.vtu")
        
    print("")
    print("--- %s seconds ---" % (time.time() - start_time))
    print('done')