def push(particles,dt):
	for particle_set in particles:
			for particle in particle_set:
				particle.position_old = particle.position

				particle.position[0] = particle.position[0] + particle.velocity[0]*dt
				particle.position[1] = particle.position[1] + particle.velocity[1]*dt
				particle.position[2] = particle.position[2] + particle.velocity[2]*dt