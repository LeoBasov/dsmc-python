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
	print(40*"-")
	print('Execution finished')
	print(40*"=")

def loop(dt,itters,particles):
	for i in range(itters):
		print('Itterartion {} complete'.format(i + 1))

def generate_particles(input_values):
	particles = []

	for set in input_values.plasma:
		particle_loc = generator.generate_particles(set.number,set,set.domain)

		particles.append(particle_loc)

	return particles

def main():
	print_header()
	
	input_values = read_xml(sys.argv[1])
	particles = generate_particles(input_values)

	loop(input_values.time.dt,input_values.time.itters,particles)

	print_footer()

main()