
class Cuboid:
	"""Cuboid shaped domain"""
	def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax):
		self.xmin = xmin
		self.xmax = xmax
		self.ymin = ymin
		self.ymax = ymax
		self.zmin = zmin
		self.zmax = zmax

	def check_if_inside(self, position):
		"""Function used to check if position lies inside of domain"""
		inx = (position[0] <= self.xmax) and (position[0] >= self.xmin)
		iny = (position[1] <= self.ymax) and (position[1] >= self.ymin)
		inz = (position[2] <= self.zmax) and (position[2] >= self.zmin)

		return inx and iny and inz