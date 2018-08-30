from plasma import particle
from geometry import domain
import random

def generate_particles(number,values):
	particles = []

	for i in range(number):
		part = particle.Particle()

		part.weight = values.weight
		part.mass = values.mass

		particles.append(part)

	return particles