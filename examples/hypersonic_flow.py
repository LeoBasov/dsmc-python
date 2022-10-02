import dsmc
import time
import numpy as np

if __name__ == '__main__':
    # general parameters
    solver = dsmc.DSMC()
    domain = [(-3.0, 3.0), (-1.5, 1.5), (-0.025, 0.025)]
    obj = [(-0.25, 0.25), (-0.25, 0.25), (-0.5, 0.5)]
    dt = 1e-6
    w = 2.3e+13
    mass = 6.6422e-26
    T =  273.0
    n = 2.6e+19
    u = np.array([0.0, -3043.0, 0.0])
    niter = 1
    
    # setup solver
    solver.set_domain(domain)
    solver.set_weight(w)
    solver.set_mass(mass)
    
    # set up boundary
    solver.set_bc_type("xmin", "inflow")
    solver.set_bc_type("xmax", "inflow")
    solver.set_bc_type("ymin", "inflow")
    
    solver.set_bc_type("ymax", "open")
    
    solver.set_bc_values("xmin", T, n, u)
    solver.set_bc_values("xmax", T, n, u)
    solver.set_bc_values("ymin", T, n, u)
    
    # set object
    solver.add_object(obj)
    
    # create particles
    solver.create_particles(domain, T, n)
    
    # start timing
    start_time = time.time()
    
    for it in range(niter):
        print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)
        solver.advance(dt)
        
    print("")
    print("--- %s seconds ---" % (time.time() - start_time))
    print('done')