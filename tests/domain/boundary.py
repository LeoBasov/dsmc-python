import dsmc.dsmc as ds
import dsmc.particles as prt
import dsmc.octree as oc
import dsmc.writer as wrt
import dsmc.boundary as bo
import dsmc.common as co
import time
import numpy as np

if __name__ == '__main__':
    # general parameters
    boundary = bo.Boundary()
    particles = prt.Particles()
    domain = np.array([(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)])
    dt = 0.0001
    w = 2.3e+18
    mass = 6.6422e-26
    T =  273.0
    n = 2.6e+19
    N = int(co.get_V(domain)*n/w)
    niter = 1000
    
    # trees
    tree_outer = oc.Octree()
    outer_pos = np.array([[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]])
    tree_outer.build(outer_pos)
    wrt.write_buttom_leafs(tree_outer, "outer_box.vtu")
    
    # setup solver
    boundary.domain = domain
    
    # create particles
    particles.create_particles(domain, mass, T, N)
    
    # start timing
    start_time = time.time()
    
    wrt.write_positions(particles, "pos_{}.csv".format(0))
    
    for it in range(niter):
        print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)
        velocities, positions, old_positions = ds._push(particles.Vel, particles.Pos, dt)
        velocities, positions, old_positions = boundary.boundary(velocities, positions, old_positions)
        
        particles.VelPos = (velocities, positions)
        wrt.write_positions(particles, "pos_{}.csv".format(it + 1))
    

    print("")
    print("--- %s seconds ---" % (time.time() - start_time))
    print('done')