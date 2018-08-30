import sys
sys.path.append('../../src')

import dsmc
import pusher
from plasma import generator
from IO.input import read_xml
from IO.output import write_particles_pos

def print_header():
	print(40*"=")
	print('DSMC python test')
	print(40*"-")

def print_footer():
	print(40*"-")
	print('Execution finished')
	print(40*"=")

def loop(dt,itters,particles,file_name,domain):
	for i in range(itters):
		pusher.push(particles,dt)

		for particle_set in particles:
			for particle in particle_set:
				domain.exec_mirrow_boundary(particle.position,particle.velocity)

		name_string = "{}_".format(i + 1)
		name_string = name_string + file_name

		write_particles_pos(name_string,particles)

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

	loop(input_values.time.dt,input_values.time.itters,particles,input_values.file_name,input_values.domain)

	print_footer()

main()