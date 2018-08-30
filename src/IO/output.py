def write_particles_pos(file_name,particles):
	with open(file_name,'w+') as file:
		file.write('#pos_x,#pos_y,#pos_z\n')

		for particle_set in particles:
			for particle in particle_set:
				file.write("{},{},{}\n".format(particle.position[0],particle.position[1],particle.position[2]))