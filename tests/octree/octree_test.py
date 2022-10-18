import dsmc.writer as wrt
import dsmc.octree as oc
import numpy as np
import argparse

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
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('N', type=int)

    args = parser.parse_args()
    
    radius = 1.0
    positions = create_particles(args.N, radius)
    octree = oc.Octree()
    octree.build(positions)
    wrt.write_buttom_leafs(octree)
    
    for i in range(len(octree.leafs)):
        box = octree.cell_boxes[i]
        leaf = octree.leafs[i]
        N = 0
        for j in range(leaf.elem_offset, leaf.elem_offset + leaf.number_elements):
            p = octree.permutations[j]
            pos = positions[p]
            
            if oc._is_inside(pos, box):
                N += 1
            else:
                raise Exception(pos, box)
                
        if N != leaf.number_elements:
            raise Exception(N, leaf.number_elements, box)
	
    with open("particles.csv", "w") as file:
        file.write("x, y, z\n")
        
        for i in range(len(octree.leafs)):
            leaf = octree.leafs[i]
            if leaf.number_children == 0 and leaf.number_elements > 0:
                for j in range(leaf.elem_offset, leaf.elem_offset + leaf.number_elements):
                    p = octree.permutations[j]
                    pos = positions[p]
                    file.write("{}, {}, {}\n".format(pos[0], pos[1], pos[2]))
	
    print("done")
