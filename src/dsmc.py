from plasma.particle import Particle
from random import shuffle
from random import uniform
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

	def _calc_collision_probability(self, number, cross_section, dt, volume, number_pairs, rel_vel, weight):
		return rel_vel*0.5*number*number*cross_section*weight*dt/(number_pairs*volume)

	def _calc_and_set_new_vels(self, particle1, particle2, rel_vel):
		v_cm = [0,0,0]
		g = [0,0,0]

		cos_chi = uniform(0.0,2.0) - 1
		sin_chi = math.sqrt(1.0 - cos_chi*cos_chi)
		epsilon = uniform(0.0,2.0)*math.pi

		v_cm[0] = (particle1.mass*particle1.velocity[0] + particle2.mass*particle2.velocity[0])/(particle1.mass + particle2.mass)
		v_cm[1] = (particle1.mass*particle1.velocity[1] + particle2.mass*particle2.velocity[1])/(particle1.mass + particle2.mass)
		v_cm[2] = (particle1.mass*particle1.velocity[2] + particle2.mass*particle2.velocity[2])/(particle1.mass + particle2.mass)

		g[0] = rel_vel*cos_chi
		g[1] = rel_vel*sin_chi*math.cos(epsilon)
		g[2] = rel_vel*sin_chi*math.sin(epsilon)

		m_1 = particle1.mass/(particle1.mass + particle2.mass)
		m_2 = particle2.mass/(particle1.mass + particle2.mass)

		particle1.velocity[0] = v_cm[0] + m_2*g[0]
		particle1.velocity[1] = v_cm[1] + m_2*g[1]
		particle1.velocity[2] = v_cm[2] + m_2*g[2]

		particle2.velocity[0] = v_cm[0] - m_1*g[0]
		particle2.velocity[1] = v_cm[1] - m_1*g[1]
		particle2.velocity[2] = v_cm[2] - m_1*g[2]

	def interact(self, particles, volume, dt):
		shuffle(particles)

		i = 1

		while i < len(particles):
			reduced_mass = self._calc_reduced_mass(particles[i - 1].mass, particles[i].mass)
			relative_velocity = self._calc_rel_vel(particles[i - 1].velocity, particles[i].velocity)
			relative_energy = 0.5*reduced_mass*relative_velocity*relative_velocity

			cross_section = self._get_cross_section(relative_energy)
			collision_probability = self._calc_collision_probability(float(len(particles)), cross_section, dt, volume, 0.5*float(len(particles)), relative_velocity, particles[i].weight)

			if collision_probability >= uniform(0.0,1.0):
				self._calc_and_set_new_vels(particles[i - 1],particles[i], relative_velocity)

			i += 2