import dsmc

if __name__ == '__main__':
    # general parameters
    solver = dsmc.DSMC()
    domain = ((-0.1e-3, 0.1e-3), (-0.1e-3, 0.1e-3), (0, 50e-3))
    dt = 1e-7
    w = 2.4134e+15
    mass = 6.6422e-26
    niter = 300
    
    # low denisty particles
    nhigh = 2.5e+20
    Thigh = 300
    Boxhigh = ((-0.1e-3, 0.1e-3), (-0.1e-3, 0.1e-3), (25e-3, 50e-3))
    
    # high denisty particles
    nlow = 2.5e+19
    Tlow = 300
    Boxlow = ((-0.1e-3, 0.1e-3), (-0.1e-3, 0.1e-3), (0, 25e-3))
    
    solver.set_domain(domain)
    solver.set_weight(w)
    solver.set_mass(mass)
    
    solver.create_particles(Boxlow, Tlow, nlow)
    solver.create_particles(Boxhigh, Thigh, nhigh)
    
    with open("test.csv", "w") as file:
        for it in range(niter):
            print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)
            solver.advance(dt)
            
            for pos in solver.particles.Pos:
                file.write("{:.4e},".format(pos[2]))
                
            file.write("\n")

    print("")
    print('done')
