import numpy as np
import dsmc.diagnostics as dia
import time

def create_particles(N, radius):
	positions = np.empty((N + 1, 3))

	for i in range(N):
		phi = np.random.random() * 2.0 * np.pi
		theta = np.random.random() * np.pi
		r = np.random.normal(0.0, 0.01)
		theta1 = np.random.random() * np.pi - 0.5 * np.pi
		dis = np.array((r * np.sin(theta) * np.cos(phi), r * np.sin(theta) * np.sin(phi), r * np.cos(theta)))
		offset = np.array((radius * np.sin(theta1), radius * np.cos(theta1), 0.0))

		dis += offset;
		positions[i] = dis

	positions[N] = np.array((0.0, -1.0, 0.0))
	
	return positions

if __name__ == "__main__":
    N = 100000
    radius = 1.0
    positions = create_particles(N, radius)
    
    # set up sort
    Nx = 10
    Ny = 10
    x0 = -1.5
    x1 = 1.5
    y0 = -1.5
    y1 = 1.5
    
    start_time = time.time()
    permutations, leafs, boxes = dia._sort_2d(positions, x0, x1, y0, y1, Nx, Ny)
    print("--- %s seconds ---" % (time.time() - start_time))   
    print("done")