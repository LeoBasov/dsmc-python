import dsmc

if __name__ == '__main__':
    # general parameters
    solver = dsmc.DSMC()
    domain = ((-0.1e-3, 0.1e-3), (-0.1e-3, 0.1e-3), (0, 50e-3))
    dt = 1e-7
    w = 2.4134e+10
    mass = 6.6422e-26
    niter = 300
    
    # low denisty particles
    nhigh = 2.5e+14
    Thigh = 300
    Boxhigh = ((-0.1e-3, 0.1e-3), (-0.1e-3, 0.1e-3), (25e-3, 50e-3))
    
    # high denisty particles
    nlow = 2.5e+13
    Tlow = 300
    Boxlow = ((-0.1e-3, 0.1e-3), (-0.1e-3, 0.1e-3), (0, 25e-3))
    
    solver.set_domain(domain)
    solver.set_weight(w)
    
    solver.create_particles(Boxlow, mass, Tlow, nlow)
    solver.create_particles(Boxhigh, mass, Thigh, nhigh)
    
    for it in range(niter):
        print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)

    print("")
    print('done')
