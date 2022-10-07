import dsmc.dsmc as ds
import dsmc.particles as prt
import dsmc.octree as oc
import dsmc.writer as wrt
import dsmc.boundary as bo
import dsmc.common as co
import time
import numpy as np

def create_pos_and_vels():
    positions = np.zeros((6, 3))
    velocitiies = np.zeros((6, 3))
    
    # x
    positions[0][0] = -0.75
    velocitiies[0][0] = 100.0
    
    positions[1][0] = 0.75
    velocitiies[1][0] = -100.0
    
    # y
    positions[2][1] = -0.75
    velocitiies[2][1] = 100.0
    
    positions[3][1] = 0.75
    velocitiies[3][1] = -100.0
    
    # z
    positions[4][2] = -0.75
    velocitiies[4][2] = 100.0
    
    positions[5][2] = 0.75
    velocitiies[5][2] = -100.0
    
    return (velocitiies, positions)

if __name__ == '__main__':
    # general parameters
    boundary = bo.Boundary()
    particles = prt.Particles()
    domain = np.array([(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)])
    dt = 0.0001
    w = 2.3e+16
    mass = 6.6422e-26
    T =  273.0
    n = 2.6e+19
    N = int(co.get_V(domain)*n/w)
    niter = 1000
    
    # trees
    tree_outer = oc.Octree()
    outer_pos = np.array([[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]])
    tree_outer.build(outer_pos)
    #wrt.write_buttom_leafs(tree_outer, "outer_box.vtu")
    
    # setup solver
    boundary.domain = domain
    
    # create particles
    particles.create_particles(domain, mass, T, N)
    
    velocities, positions = particles.VelPos
        
    #particles.VelPos = create_pos_and_vels()
    
    # start timing
    start_time = time.time()
    
    wrt.write_positions(particles, "pos_{}.csv".format(0))
    
    E0 = 0.0
    
    for vel in particles.Vel:
        E0 += vel.dot(vel)
    
    for it in range(niter):
        E = 0.0
        
        velocities, positions, old_positions = ds._push(particles.Vel, particles.Pos, dt)
        velocities, positions, old_positions = boundary.boundary(velocities, positions, old_positions)
        
        particles.VelPos = (velocities, positions)
        wrt.write_positions(particles, "pos_{}.csv".format(it + 1))
        
        for vel in particles.Vel:
            E += vel.dot(vel)
        
        print("iteration {:4}/{}, N particles {}/{}, Efrac {}".format(it + 1, niter, particles.N, N, E/E0), end="\r", flush=True)
    

    print("")
    print("--- %s seconds ---" % (time.time() - start_time))
    print('done')