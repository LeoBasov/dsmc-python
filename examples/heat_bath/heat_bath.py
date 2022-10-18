import dsmc
import dsmc.particles as prt
import matplotlib.pyplot as plt
import numpy as np
import math
import time

def maxwell(x, T):
    return 2.0 * np.sqrt(x) * np.exp(-x/T) / (math.pow(T, 3.0/2.0) * np.sqrt(math.pi))
    
def calc_x(velocities, mass):
    return np.array([mass*vel.dot(vel)/(2.0*prt.kb) for vel in velocities])

if __name__ == '__main__':
    # general parameters
    solver = dsmc.DSMC()
    domain = ((-1.0e-3, 1.0e-3), (-1.0e-3, 1.0e-3), (-1.0e-3, 1.0e-3))
    dt = 1e-5
    w = 0.5e+8
    mass = 6.6422e-26
    niter = 100
    Nbins = 200
    
    # particles
    n = 1.0e+20
    T = 300
    u = np.array((1000.0, 0.0, 0.0))
    
    
    solver.set_domain(domain)
    solver.set_weight(w)
    solver.set_mass(mass)
    
    solver.create_particles(domain, T, n, u)
    
    # time
    start_time = time.time()
    
    # set up plot
    xmax = 15000
    
    Tnew = prt.calc_temperature(solver.particles.Vel, solver.mass)
    xm = np.linspace(0, xmax, 1000)
    dist = [maxwell(xi, Tnew) for xi in xm]
    x = calc_x(solver.particles.Vel, solver.mass)
    
    fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(16, 6))
    
    ax0.plot(xm, dist)
    ax0.hist(x, Nbins, density=True)
    ax0.set_ylabel("probabilty density")
    ax0.set_xlabel("m v^2 / (2 kb)")
    ax0.set_title("initial condition")
    ax0.set_xlim([0, xmax])
    ax0.set_ylim([0, 0.00040])
    
    for it in range(niter):
        print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)
        solver.advance(dt)
        x = calc_x(solver.particles.Vel, solver.mass)
    
    print("")
    print("--- %s seconds ---" % (time.time() - start_time))
    
    ax1.plot(xm, dist)
    ax1.hist(x, Nbins, density=True)
    ax1.set_ylabel("probabilty density")
    ax1.set_xlabel("m v^2 / (2 kb)")
    ax1.set_title("final condition")
    ax1.set_xlim([0, xmax])
    ax1.set_ylim([0, 0.00040])
    
    fig.suptitle("Argon Heat Bath")
    plt.show() 
    
    print('done')
