"""Module designed to read files and present them in internaly readable format"""

import xml.etree.ElementTree as ET


def read_domain(node):
	pass

def read_xml(file_name):
	try:
		tree = ET.parse(file_name)
		root  = tree.getroot()

		for child in root:
			print(child.tag, child.attrib)
	except:
		raise
	else:
		pass
	finally:
		pass