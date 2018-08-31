from plasma.particle import Particle
from random import shuffle

class DSMC:
	def __init__(self, cross_sections):
		self.cross_sections = cross_sections

	def interact(self, particles):
		shuffle(particles)