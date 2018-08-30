from plasma import particle
from geometry import domain
import random

class Values:
	def __init__(self,weight,mass,temperature):
		self.weight = weight
		self.mass = mass
		self.temperature = temperature

def generate_particles(number,values,domain):
	particles = []

	for i in range(number):
		part = particle.Particle()

		part.weight = values.weight
		part.mass = values.mass

		part.position[0] = random.uniform(domain.xmin,domain.xmax)
		part.position[1] = random.uniform(domain.ymin,domain.ymax)
		part.position[2] = random.uniform(domain.zmin,domain.zmax)

		particles.append(part)

	return particles