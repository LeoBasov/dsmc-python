import dsmc
import dsmc.diagnostics as dia

def write2file(solver, n_file, T_file, p_file):
    bins, box, x = dia.sort_bin(solver.particles.Pos, 2, Nbins)
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

if __name__ == '__main__':
    # general parameters
    solver = dsmc.DSMC()
    domain = [(-0.0001, 0.0001), (-0.0001, 0.0001), (0.0, 0.1)]
    dt = 1e-7
    w = 1e+9
    mass = 6.6422e-26
    niter = 300
    Nbins = 100
    
    # high denisty particles
    nhigh = 2.41432e+22 
    Thigh = 300
    Boxhigh = [(-0.0001, 0.0001), (-0.0001, 0.0001), (0.0, 0.05)]
    
    # low denisty particles
    nlow = 2.41432e+21
    Tlow = 300
    Boxlow = [(-0.0001, 0.0001), (-0.0001, 0.0001), (0.05, 0.1)]
    
    # setup solver
    solver.set_domain(domain)
    solver.set_weight(w)
    solver.set_mass(mass)
    
    solver.create_particles(Boxlow, Tlow, nlow)
    solver.create_particles(Boxhigh, Thigh, nhigh)
    
    # open files
    n_file = open("n.csv", "w")
    T_file = open("T.csv", "w")
    p_file = open("p.csv", "w")
    
    for it in range(niter):
        print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)
        solver.advance(dt)        
        write2file(solver, n_file, T_file, p_file)

	# close files
    n_file.close()
    T_file.close()
    p_file.close()        
    
    print("")
    print('done')
