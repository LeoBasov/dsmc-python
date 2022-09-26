import dsmc
import dsmc.writer as wrt

if __name__ == '__main__':
	# general parameters
    solver = dsmc.DSMC()
    domain = [(-0.5, 0.5), (-0.5, 0.5), (-0.5, 0.5)]
    dt = 1e-5
    w = 1e+17
    n = 1e18
    mass = 6.6422e-26
    T = 300
    niter = 300
    
    # setup solver
    solver.set_domain(domain)
    solver.set_weight(w)
    solver.set_mass(mass)
    
    solver.create_particles(domain, T, n)
    
    for it in range(niter):
        print("iteration {:4}/{:4}".format(it + 1, niter), end="\r", flush=True)
        solver.advance(dt, collisions=False)
        
        wrt.write_buttom_leafs(solver.octree, "octree_{}.vtu".format(it))
        
        with open("particles_{}.csv".format(it), "w") as file:
            file.write("x, y, z\n")
            for pos in solver.particles.Pos:
                file.write("{},{},{}\n".format(pos[0], pos[1], pos[2]))
            
    
    print("done")