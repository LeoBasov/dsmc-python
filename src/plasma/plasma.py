"""
Plasma Module
Author: Leo Basov

This module conttains all classes and algorithms associated with the representation of plasma
"""

class Particle:
	"""Main particle class"""
	def __init__(self):
		self.mass = 1.0
		self.weight = 1.0
		self.species = '40Ar'
		self.velocity = [0.0,0.0,0.0]
		self.position = [0.0,0.0,0.0]