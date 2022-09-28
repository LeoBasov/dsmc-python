import dsmc
import dsmc.octree as oc
import dsmc.writer as wrt
import numpy as np

if __name__ == '__main__':
	# general parameters
    solver = dsmc.DSMC()
    tree = oc.Octree()
    domain = [(-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5)]
    positions = np.array([(-0.5, -0.5, -0.5), (0.5, 0.5, 0.5)])
    dt = 1e-5
    w = 1e+16
    n = 1e18
    mass = 6.6422e-26
    T = 300
    niter = 1000
    
    # setup solver
    solver.set_domain(domain)
    solver.set_weight(w)
    solver.set_mass(mass)
    
    solver.create_particles(domain, T, n)
    
    solver.set_boundary("xmax", "open")
    
    solver.set_boundary("ymax", "open")
    solver.set_boundary("ymin", "open")
    solver.set_boundary("zmax", "open")
    solver.set_boundary("zmin", "open")
    
    solver.set_boundary("xmin", "inflow")
    
    solver.boundary.u[0][0] = np.array([1000, 0, 0])
    
    tree.build(positions)
    
    for it in range(niter):
        print("iteration {:4}/{:4}".format(it + 1, niter), end="\r", flush=True)
        
        try:
            solver.advance(dt, collisions=False, octree=False)
        except Exception as e:
            print(e)
            break
        
        with open("particles_{}.csv".format(it), "w") as file:
            file.write("x, y, z\n")
            for pos in solver.particles.Pos:
                file.write("{},{},{}\n".format(pos[0], pos[1], pos[2]))
                
    wrt.write_buttom_leafs(tree, "octree.vtu")
            
    
    print("done")