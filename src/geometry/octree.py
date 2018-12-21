"""Octree class for pair selection for dsmc."""

from geometry.domain import Cuboid
import math

class Leaf:
	"""Info"""
	def __init__(self, domain):
		self.domain = domain
		self.particles = []
		self.maximum_cross_section = 3.30606633873877e-19

		self.parent = None
		self.children = 8*[None]

	def __iter__(self):
		return self.children.__iter__()

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

	def _calc_mean_free_path(self):
		number_dens = self._calc_number_dens()

		if number_dens > 0.0:
			return 1.0/(number_dens*self.maximum_cross_section)
		else:
			return math.inf

	def check_resultion_criterion(self):
		if len(self.particles) > 1:
			return True
		else:
			return False

	def has_children(self):
		return not (self.children[0] == None)



class Octree:
	"""docstring for Octree"""
	def __init__(self):
		self.leafs = []
		self.buttom_leafs = []

	def _subdivide(self, domain, parent):
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

		parent.children = new_leafs

		return new_leafs

	def _build_next_level(self, leaf):
		loc_leafs = self._subdivide(leaf.domain, leaf)
		particles = leaf.particles

		for loc_leaf in loc_leafs:
			self.leafs.append(loc_leaf)

			particles = loc_leaf.sort(particles)

			if loc_leaf.check_resultion_criterion():
				self._build_next_level(loc_leaf)
			else:
				self.buttom_leafs.append(loc_leaf)

	def build(self, particles, domain):
		self.leafs = []
		self.buttom_leafs = []

		main_leaf = Leaf(domain)
		main_leaf.sort(particles)

		self.leafs.append(main_leaf)

		if main_leaf.check_resultion_criterion():
			self._build_next_level(main_leaf)
		else:
			self.buttom_leafs.append(main_leaf)
