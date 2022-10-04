import dsmc
import dsmc.octree as oc
import dsmc.writer as wrt
import time
import numpy as np

def create_pos_and_vels():
    positions = np.zeros((6, 3))
    velocitiies = np.zeros((6, 3))
    
    # x
    positions[0][0] = -0.75
    velocitiies[0][0] = 1.0
    
    positions[1][0] = 0.75
    velocitiies[1][0] = -1.0
    
    # y
    positions[2][1] = -0.75
    velocitiies[2][1] = 1.0
    
    positions[3][1] = 0.75
    velocitiies[3][1] = -1.0
    
    # z
    positions[4][2] = -0.75
    velocitiies[4][2] = 1.0
    
    positions[5][2] = 0.75
    velocitiies[5][2] = -1.0
    
    return (velocitiies, positions)

def write_partices(positions, it):
    with open("particles_{}.csv".format(it), "w") as file:
        file.write("x, y, z\n")
        for pos in positions:
            file.write("{}, {}, {}\n".format(pos[0], pos[1], pos[2]))

if __name__ == '__main__':
    # general parameters
    solver = dsmc.DSMC()
    domain = [(-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0)]
    obj = [(-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5)]
    dt = 0.01
    w = 2.3e+14
    mass = 6.6422e-26
    T =  273.0
    n = 2.6e+19
    niter = 200
    
    # trees
    tree_inner = oc.Octree()
    tree_outer = oc.Octree()
    
    inner_pos = np.array([[-0.5, -0.5, -0.5], [0.5, 0.5, 0.5]])
    outer_pos = np.array([[-1.0, -1.0, -1.0], [1.0, 1.0, 1.0]])
    
    tree_inner.build(inner_pos)
    tree_outer.build(outer_pos)
    
    # setup solver
    solver.set_domain(domain)
    solver.set_weight(w)
    solver.set_mass(mass)
    
    # set object
    solver.add_object(obj)
    
    # create particles
    solver.particles.VelPos = create_pos_and_vels()
    
    # start timing
    start_time = time.time()
    
    for it in range(niter):
        print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)
        solver.advance(dt, octree=False)
        write_partices(solver.particles.Pos, it)
        
    wrt.write_buttom_leafs(tree_inner, "inner_box.vtu")
    wrt.write_buttom_leafs(tree_outer, "outer_box.vtu")

    print("")
    print("--- %s seconds ---" % (time.time() - start_time))
    print('done')