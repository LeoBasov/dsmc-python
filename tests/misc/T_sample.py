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
	solver = dsmc.DSMC()
	domain = ((-0.1e-3, 0.1e-3), (-0.1e-3, 0.1e-3), (0, 50e-3))
	w = 1e6
	mass = 6.6422e-26
	T = 300
	n = 1e+20
	Nbins = 100
    
	solver.set_domain(domain)
	solver.set_weight(w)
	solver.set_mass(mass)
    
	solver.create_particles(domain, T, n)
    
	print(prt.calc_temperature(solver.particles.Vel, mass))
	
	x = calc_x(solver.particles.Vel, solver.mass)

	xm = np.linspace(0, 3500, 1000)
	dist = [maxwell(xi, 300) for xi in xm]
    
	plt.plot(xm, dist)
	plt.hist(x, Nbins, density=True)
	plt.show()
