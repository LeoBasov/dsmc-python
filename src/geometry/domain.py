import math


class Cuboid:
	"""Cuboid shaped domain"""
	def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax):
		self.xmin = xmin
		self.xmax = xmax
		self.ymin = ymin
		self.ymax = ymax
		self.zmin = zmin
		self.zmax = zmax

		self.volume = (self.xmax - self.xmin)*(self.ymax - self.ymin)*(self.zmax - self.zmin)
		self.diagonal = math.sqrt((self.xmax - self.xmin)**2 + (self.ymax - self.ymin)**2 + (self.zmax - self.zmin)**2)

	def check_if_inside(self, position):
		"""Function used to check if position lies inside of domain"""
		inx = (position[0] <= self.xmax) and (position[0] >= self.xmin)
		iny = (position[1] <= self.ymax) and (position[1] >= self.ymin)
		inz = (position[2] <= self.zmax) and (position[2] >= self.zmin)

		return inx and iny and inz

	def calc_factor(self, value, min, max, length, value_bigger_then_max):
		if value_bigger_then_max:
			return math.ceil((value - max)/length)
		else:
			return math.ceil((min - value)/length)

	def calc_new_value(self, value, min, max,length, value_bigger_then_max):
		if value_bigger_then_max:
			return value - self.calc_factor(value,min,max,length,value_bigger_then_max)*length
		else:
			return value + self.calc_factor(value,min,max,length,value_bigger_then_max)*length

	def calc_new_x(self, position, min, max):
		if position[0] > max:
			position[0] = self.calc_new_value(position[0],min,max,max - min,True)
		elif position[0] < min:
			position[0] = self.calc_new_value(position[0],min,max,max - min,False)

	def calc_new_y(self, position, min, max):
		if position[1] > max:
			position[1] = self.calc_new_value(position[1],min,max,max - min,True)
		elif position[1] < min:
			position[1] = self.calc_new_value(position[1],min,max,max - min,False)

	def calc_new_z(self, position, min, max):
		if position[2] > max:
			position[2] = self.calc_new_value(position[2],min,max,max - min,True)
		elif position[2] < min:
			position[2] = self.calc_new_value(position[2],min,max,max - min,False)

	def reposition(self,position):
		self.calc_new_x(position,self.xmin,self.xmax)
		self.calc_new_y(position,self.ymin,self.ymax)
		self.calc_new_z(position,self.zmin,self.zmax)

	def exec_periodic_boundary(self,position):
		if not self.check_if_inside(position):
			self.reposition(position)

	def exec_mirrow_boundary(self, position_old, position, velocity):
		while not self.check_if_inside(position):
			intersection_point = self._find_intersection_point(position_old, position)

			if position[0] < self.xmin:
				position[0] = 2*self.xmin - position[0]
				velocity[0] = (-1)*velocity[0]

			if position[1] < self.ymin:
				position[1] = 2*self.ymin - position[1]
				velocity[1] = (-1)*velocity[1]

			if position[2] < self.zmin:
				position[2] = 2*self.zmin - position[2]
				velocity[2] = (-1)*velocity[2]

			if position[0] > self.xmax:
				position[0] = 2*self.xmax - position[0]
				velocity[0] = (-1)*velocity[0]

			if position[1] > self.ymax:
				position[1] = 2*self.ymax - position[1]
				velocity[1] = (-1)*velocity[1]

			if position[2] > self.zmax:
				position[2] = 2*self.zmax - position[2]
				velocity[2] = (-1)*velocity[2]

	def _find_intersection_point(self, position_old, position):
		if (self.check_if_inside(position_old) != self.check_if_inside(position)):
			line_factor = [None]*6
			intersection_point = [None]*3

			line_factor[0] = (self.xmin - position_old[0])/(position[0] - position_old[0])
			line_factor[1] = (self.xmax - position_old[0])/(position[0] - position_old[0])

			line_factor[2] = (self.ymin - position_old[1])/(position[1] - position_old[1])
			line_factor[3] = (self.ymax - position_old[1])/(position[1] - position_old[1])

			line_factor[4] = (self.zmin - position_old[2])/(position[2] - position_old[2])
			line_factor[5] = (self.zmax - position_old[2])/(position[2] - position_old[2])

			for i in range(len(line_factor)):
				line_factor[i] = abs(line_factor[i])

			factor = min(line_factor)

			intersection_point[0] = position_old[0] + factor*(position[0] - position_old[0])
			intersection_point[1] = position_old[1] + factor*(position[1] - position_old[1])
			intersection_point[2] = position_old[2] + factor*(position[2] - position_old[2])

			return intersection_point
		else:
			raise Exception('Both positions are either in or out',position_old,position)