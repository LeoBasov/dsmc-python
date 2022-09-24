import dsmc
import dsmc.particles as prt
import matplotlib.pyplot as plt
import numpy as np
import math

def maxwell(x, T):
    return 2.0 * np.sqrt(x) * np.exp(-x/T) / (math.pow(T, 3.0/2.0) * np.sqrt(math.pi))
    
def calc_x(velocities, mass):
    return np.array([mass*vel.dot(vel)/(2.0*prt.kb) for vel in velocities])

if __name__ == '__main__':
    # general parameters
    solver = dsmc.DSMC()
    domain = ((-0.1e-3, 0.1e-3), (-0.1e-3, 0.1e-3), (0, 50e-3))
    dt = 1e-5
    w = 2.4134e+7
    mass = 6.6422e-26
    niter = 500
    Nbins = 100
    
    # particles
    n = 2.5e+20
    T = 300
    Box = ((-0.1e-3, 0.1e-3), (-0.1e-3, 0.1e-3), (25e-3, 50e-3))
    u = np.array((1000.0, 0.0, 0.0))
    
    
    solver.set_domain(domain)
    solver.set_weight(w)
    solver.set_mass(mass)
    
    solver.create_particles(Box, T, n, u)
    
    for it in range(niter):
        print("iteration {:4}/{}".format(it + 1, niter), end="\r", flush=True)
        solver.advance(dt)
        x = calc_x(solver.particles.Vel, solver.mass)

    xm = np.linspace(0, 3500, 1000)
    dist = [maxwell(xi, 300) for xi in xm]
    
    plt.plot(xm, dist)
    plt.hist(x, Nbins, density=True)
    plt.show()  
    
    print("")
    print('done')
