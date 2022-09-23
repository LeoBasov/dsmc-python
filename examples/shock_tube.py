import dsmc
import dsmc.diagnostics as dia
import matplotlib.pyplot as plt

if __name__ == '__main__':
    # general parameters
    solver = dsmc.DSMC()
    domain = ((-0.1e-3, 0.1e-3), (-0.1e-3, 0.1e-3), (0, 50e-3))
    dt = 1e-7
    w = 2.4134e+7
    mass = 6.6422e-26
    niter = 300
    Nbins = 100
    
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
    
    n_file = open("n.csv", "w")
    T_file = open("T.csv", "w")
    p_file = open("p.csv", "w")
    
    for it in range(niter):
        print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)
        solver.advance(dt)
            
        bins, box = dia.sort_bin(solver.particles.Pos, 2, Nbins)
        N = [len(b) for b in bins]
        n = dia.calc_n(bins, box, 2, solver.w)
        T = dia.calc_T(bins, solver.particles.Vel, mass)
        p = dia.calc_p(n, T)
      
        for i in range(Nbins):
            n_file.write("{},".format(n[i]))
            T_file.write("{},".format(T[i]))
            p_file.write("{},".format(p[i]))
                
        n_file.write("\n")
        T_file.write("\n")
        p_file.write("\n")

    n_file.close()
    T_file.close()
    p_file.close()        
    
    print("")
    print('done')
