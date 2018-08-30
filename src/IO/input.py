"""Module designed to read files and present them in internaly readable format"""

import xml.etree.ElementTree as ET
from geometry.domain import Cuboid

class Empty:
	pass


def read_domain(node):
	geometry = node.find("geometry")

	x_min = float(geometry.get('x_min'))
	x_max = float(geometry.get('x_max'))

	y_min = float(geometry.get('y_min'))
	y_max = float(geometry.get('y_max'))

	z_min = float(geometry.get('z_min'))
	z_max = float(geometry.get('z_max'))

	return Cuboid(x_min, x_max, y_min, y_max, z_min, z_max)

def read_plasma(node):
	retValls = []
	sets = node.findall('set')

	for set in sets:
		val = Empty()

		weight_node = set.find("weight")
		mass_node = set.find("mass")
		temperatur_node = set.find("temperature")
		drift_velocity_node = set.find("drift_velocity")
		domain_node = set.find("domain")

		val.weight = float(weight_node.get('value'))
		val.mass = float(mass_node.get('value'))
		val.temperature = [float(temperatur_node.get('x_coord')),float(temperatur_node.get('y_coord')),float(temperatur_node.get('z_coord'))]
		val.drift_velocity = [float(drift_velocity_node.get('x_coord')),float(drift_velocity_node.get('y_coord')),float(drift_velocity_node.get('z_coord'))]
		val.domain = read_domain(domain_node)

		retValls.append(val)

	return retValls


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

	except:
		raise
	else:
		pass
	finally:
		pass