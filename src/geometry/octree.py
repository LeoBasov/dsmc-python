"""Octree class for pair selection for dsmc."""

from geometry.domain import Cuboid

class Leaf:
	"""Info"""
	def __init__(self, domain):
		self.domain = domain
		self.particles = []

	def sort(self,particles):
		retParticles = []

		for particle in particles:
			if self.domain.check_if_inside(particle.position):
				self.particles.append(particle)
			else:
				retParticles.append(particle)

		return retParticles

	def _calc_number_dens(self):
		tot_particle_number = 0.0

		for particle in self.particles:
			tot_particle_number += particle.weight

		return tot_particle_number/self.domain.volume



class Octree:
	"""docstring for Octree"""
	def __init__(self):
		self.leafs = []
		self.buttom_leafs = []

	def build(self, particles, domain):
		pass

	def _subdivide(self, domain):
		new_leafs = []

		x_half = 0.5*(domain.xmin + domain.xmax)
		y_half = 0.5*(domain.ymin + domain.ymax)
		z_half = 0.5*(domain.zmin + domain.zmax)

		child_geo1 = Cuboid(x_half, domain.xmax, y_half,domain.ymax, z_half, domain.zmax)
		child_geo2 = Cuboid(domain.xmin, x_half, y_half,domain.ymax, z_half, domain.zmax)
		child_geo3 = Cuboid(domain.xmin, x_half, domain.ymin,y_half, z_half, domain.zmax)
		child_geo4 = Cuboid(x_half, domain.xmax, domain.ymin,y_half, z_half, domain.zmax)

		child_geo5 = Cuboid(x_half, domain.xmax, y_half,domain.ymax, domain.zmin, z_half)
		child_geo6 = Cuboid(domain.xmin, x_half, y_half,domain.ymax, domain.zmin, z_half)
		child_geo7 = Cuboid(domain.xmin, x_half, domain.ymin,y_half, domain.zmin, z_half)
		child_geo8 = Cuboid(x_half, domain.xmax, domain.ymin,y_half, domain.zmin, z_half)

		new_leafs.append(Leaf(child_geo1))
		new_leafs.append(Leaf(child_geo2))
		new_leafs.append(Leaf(child_geo3))
		new_leafs.append(Leaf(child_geo4))

		new_leafs.append(Leaf(child_geo5))
		new_leafs.append(Leaf(child_geo6))
		new_leafs.append(Leaf(child_geo7))
		new_leafs.append(Leaf(child_geo8))

		return new_leafs
