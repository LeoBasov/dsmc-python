import sys
sys.path.append('../../src')

from dsmc import DSMC
import pusher
from plasma import generator
from IO.input import read_xml
from IO.output import write_particles_pos
from geometry.octree import Octree

def diagnose(number_bins, z_min, z_max, particles,area):
	dict_part = {}
	bin_size = (z_max - z_min)/number_bins
	pos = z_min
	weight = particles[0].weight
	volume = area*bin_size

	denisties = []

	for particle in particles:
		for particle in particle_set:
			dict_part[particle.position[2]] = particle

	sorted_list = sorted(list(dict_part))
	added = 0

	for part_pos in sorted_list:
		if part_pos < pos + bin_size:
			added += weight
		else:
			denisties.append(added/volume)
			added = 0
			pos += bin_size

	return denisties

def print_header():
	print(40*"=")
	print('DSMC python test')
	print(40*"-")

def print_footer():
	print(40*"-")
	print('Execution finished')
	print(40*"=")

def print_praView_file(file_name,iter,particles):
	name_string = "{}_".format(iter + 1)
	name_string = name_string + file_name

	write_particles_pos(name_string,particles)

def execute_mirrow_boundary(particles,domain):
	for particle_set in particles:
		for particle in particle_set:
			domain.exec_mirrow_boundary(particle.position,particle.velocity)

def create_leafs(particles,domain):
	particles_loc = []
	tree = Octree()

	for particle_set in particles:
		for particle in particle_set:
			particles_loc.append(particle)

	tree.build(particles_loc,domain)

	return tree.buttom_leafs


def loop(dt,itters,particles,file_name,domain,cross_sections):
	dsmc = DSMC(cross_sections)

	for i in range(itters):
		leafs = create_leafs(particles,domain)

		for leaf in leafs:
			dsmc.interact(leaf.particles, leaf.domain.volume, dt)

		pusher.push(particles,dt)
		execute_mirrow_boundary(particles,domain)

		print_praView_file(file_name,i,particles)
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

	loop(input_values.time.dt,input_values.time.itters,particles,input_values.file_name,input_values.domain,input_values.cross_sections)

	print_footer()

main()