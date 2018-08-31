from plasma.particle import Particle
from random import shuffle
import math

class DSMC:
	def __init__(self, cross_sections):
		self.cross_sections = cross_sections
		self.sorted_keys= sorted(list(self.cross_sections))

	def _calc_reduced_mass(self, mass1, mass2):
		return mass1*mass2/(mass1 + mass2)

	def _calc_rel_vel(self, vel1, vel2):
		g = [0,0,0]

		g[0] = vel1[0] - vel2[0]
		g[1] = vel1[1] - vel2[1]
		g[2] = vel1[2] - vel2[2]

		return math.sqrt(g[0]*g[0] + g[1]*g[1] + g[2]*g[2])

	def _interpolate(self, min, max, val):
		return self.cross_sections[min]*(1.0 - ((val - min)/(max - min))) + self.cross_sections[max]*((val - min)/(max - min))

	def _get_cross_section(self, energy):
		cross_section = 0.0
		i = 1

		while i < len(self.sorted_keys):
			if energy >= self.sorted_keys[i - 1] and energy <= self.sorted_keys[i]:
				cross_section = self._interpolate(self.sorted_keys[i - 1], self.sorted_keys[i],  energy)
				break

			i += 2

		return cross_section

	def interact(self, particles, volume, dt):
		shuffle(particles)

		i = 1

		while i < len(particles):
			reduced_mass = self._calc_reduced_mass(particles[i - 1].mass, particles[i].mass)
			relative_velocity = self._calc_rel_vel(particles[i - 1].velocity, particles[i].velocity)
			relative_energy = 0.5*reduced_mass*relative_velocity*relative_velocity

			cross_sections = self._get_cross_section(relative_energy)

			i += 2