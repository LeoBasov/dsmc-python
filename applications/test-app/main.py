import sys
sys.path.append('../../src')

from dsmc import DSMC
import pusher
from plasma import generator
from IO.input import read_xml
from IO.output import write_particles_pos
from geometry.octree import Octree

def diagnose(number_bins, z_min, z_max, particles,area):
	print('Started shock_tube diagnostic')
	dict_part = {}
	bin_size = (z_max - z_min)/number_bins
	pos = z_min
	weight = particles[0][0].weight
	volume = area*bin_size

	denisties = []

	for particle_set in particles:
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
	print('Execution finished')
	print(40*"=")

def print_praView_file(file_name,iter,particles):
	name_string = "{}_".format(iter + 1)
	name_string = name_string + file_name

	write_particles_pos(name_string,particles)

def execute_boundary(particles, domain):
	if domain.type == 'specular_scattering':
		execute_mirrow_boundary(particles, domain)
	elif domain.type == 'diffuse_scattering':
		execute_diffuse_boundary(particles, domain)

def execute_mirrow_boundary(particles,domain):
	print('Started mirror boundary')
	for particle_set in particles:
		for particle in particle_set:
			domain.exec_mirrow_boundary(particle.position_old, particle.position, particle.velocity)

def execute_diffuse_boundary(particles,domain):
	print('Started diffuse boundary')
	for particle_set in particles:
		for particle in particle_set:
			domain.exec_diffuse_scattering(particle.position_old, particle.position, particle.velocity, particle.mass)

def create_leafs(particles,domain):
	particles_loc = []
	tree = Octree()

	for particle_set in particles:
		for particle in particle_set:
			particles_loc.append(particle)

	tree.build(particles_loc,domain)

	return tree


def loop(dt,itters,particles,file_name,domain,cross_sections):
	file = open("number_densities.csv",'w+')

	dsmc = DSMC(cross_sections)

	for i in range(itters):
		print('Started pusher')
		pusher.push(particles,dt)

		execute_boundary(particles,domain)

		print('Started building octree')
		tree = create_leafs(particles,domain)

		print('Started dsmc')
		first_leaf = tree.leafs[0]

		_execute_next_level(first_leaf, dt, dsmc)

		number_densities = diagnose(1000, domain.zmin, domain.zmax,particles,(domain.xmax - domain.xmin)*(domain.ymax - domain.ymin))

		print('Stared writing number deinsities')

		for density in number_densities:
			file.write('{},'.format(density))

		file.write('\n')

		#print_praView_file(file_name,i,particles)

		print('Itterartion {} complete'.format(i + 1))
		print(40*"-")

	file.close()

def _execute_next_level(leaf, dt, dsmc):
	particles = []
	volume = 0.0

	if not leaf.has_children():
		return

	for child in leaf:
		if len(child.particles) == 1:
			particles += child.particles
			volume += child.domain.volume

	dsmc.interact(particles, volume, dt)

	for child in leaf:
		if len(child.particles) > 1:
			_execute_next_level(child, dt, dsmc)

def generate_particles(input_values):
	particles = []

	for set in input_values.plasma:
		particle_loc = generator.generate_particles(set.number,set,set.domain)

		particles.append(particle_loc)

	return particles

def main():
	print_header()
	
	print("Reading input file")
	input_values = read_xml(sys.argv[1])
	print("Finished reading")
	print(40*"-")

	print('Generating particles')
	particles = generate_particles(input_values)
	print("Finished generating")
	print(40*"-")

	loop(input_values.time.dt,input_values.time.itters,particles,input_values.file_name,input_values.domain,input_values.cross_sections)

	print_footer()

main()