from plasma import particle
import random
import math

def generate_particles(number, values, domain):
	"""Function used to generate particles"""
	particles = []
	boltzmann_const = 1.38064852e-23

	for i in range(number):
		part = particle.Particle()

		part.weight = values.weight
		part.mass = values.mass

		part.velocity[0] = random.gauss(0.0, math.sqrt(boltzmann_const * values.temperature[0]/values.mass)) + values.drift_velocity[0]
		part.velocity[1] = random.gauss(0.0, math.sqrt(boltzmann_const * values.temperature[1]/values.mass)) + values.drift_velocity[1]
		part.velocity[2] = random.gauss(0.0, math.sqrt(boltzmann_const * values.temperature[2]/values.mass)) + values.drift_velocity[2]

		part.position[0] = random.uniform(domain.xmin, domain.xmax)
		part.position[1] = random.uniform(domain.ymin, domain.ymax)
		part.position[2] = random.uniform(domain.zmin, domain.zmax)

		particles.append(part)

	return particles