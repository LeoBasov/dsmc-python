"""Module designed to read files and present them in internaly readable format"""

import xml.etree.ElementTree as ET
from geometry.domain import Cuboid

class Empty:
	pass


def read_domain(node):
	geometry = node.find("geometry")
	boundary = node.find("boundary")

	x_min = float(geometry.get('x_min'))
	x_max = float(geometry.get('x_max'))

	y_min = float(geometry.get('y_min'))
	y_max = float(geometry.get('y_max'))

	z_min = float(geometry.get('z_min'))
	z_max = float(geometry.get('z_max'))

	domain = Cuboid(x_min, x_max, y_min, y_max, z_min, z_max)

	if type(boundary) != type(None):
		domain.type = boundary.get('type')
		domain.temperature = float(boundary.get('temperature'))
		domain.accommodation_factor = float(boundary.get('accommodation_factor'))

	return domain

def read_plasma(node):
	retValls = []
	sets = node.findall('set')

	for set in sets:
		val = Empty()

		number_node = set.find("number")
		weight_node = set.find("weight")
		mass_node = set.find("mass")
		temperatur_node = set.find("temperature")
		drift_velocity_node = set.find("drift_velocity")
		domain_node = set.find("domain")

		val.number = int(number_node.get('value'))
		val.weight = float(weight_node.get('value'))
		val.mass = float(mass_node.get('value'))
		val.temperature = [float(temperatur_node.get('x_coord')),float(temperatur_node.get('y_coord')),float(temperatur_node.get('z_coord'))]
		val.drift_velocity = [float(drift_velocity_node.get('x_coord')),float(drift_velocity_node.get('y_coord')),float(drift_velocity_node.get('z_coord'))]
		val.domain = read_domain(domain_node)

		retValls.append(val)

	return retValls

def read_time(node):
	val = Empty()

	sub_node = node.find("values")

	val.dt = float(sub_node.get('dt'))
	val.itters = int(sub_node.get('itters'))

	return val

def read_filemane(node):
	return node.find('file').get('name')

def read_dsmc(node):
	ret_vals = {}
	pairs = node.findall('pair')

	for pair in pairs:
		energy = float(pair.get('collision_energy'))
		cross_section = float(pair.get('cross_section'))

		ret_vals[energy] = cross_section

	return ret_vals


def read_xml(file_name):
	retVals = Empty()

	try:
		tree = ET.parse(file_name)
		root  = tree.getroot()

		for child in root:
			if child.attrib['name'] == 'domain':
				retVals.domain  = read_domain(child)
			elif child.attrib['name'] == 'plasma':
				retVals.plasma  = read_plasma(child)
			elif child.attrib['name'] == 'time':
				retVals.time = read_time(child)
			elif child.attrib['name'] == 'output':
				retVals.file_name = read_filemane(child)
			elif child.attrib['name'] == 'dsmc':
				retVals.cross_sections = read_dsmc(child)
	except:
		raise
	else:
		pass
	finally:
		pass

	return retVals