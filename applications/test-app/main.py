import sys
sys.path.append('../../src')

import dsmc
from plasma import generator
from IO.input import read_xml

def print_header():
	print(40*"=")
	print('DSMC python test')
	print(40*"-")

def print_footer():
	print(40*"=")

def main():
	print_header()
	
	inputValues = read_xml(sys.argv[1])
	particles = []

	for set in inputValues.plasma:
		particle_loc = generator.generate_particles(set.number,set,set.domain)

		particles.append(particle_loc)

	print_footer()

main()