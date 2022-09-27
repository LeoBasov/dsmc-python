import unittest
import numpy as np
from dsmc import octree as oc
from dsmc import particles as part
import csv

class TestOctree(unittest.TestCase):
    def test__find_bounding_box(self):
        positions = np.array([(0.0, 0.0, -1.0), (-2.0, -3.0, 0.0), (4.0, 5.0, 0.0)])
        box = oc._find_bounding_box(positions)
        
        self.assertEqual(-2.0, box[0][0])
        self.assertEqual(4.0, box[0][1])
        
        self.assertEqual(-3.0, box[1][0])
        self.assertEqual(5.0, box[1][1])
        
        self.assertEqual(-1.0, box[2][0])
        self.assertEqual(0.0, box[2][1])
        
    def test__calc_N_res(self):
        w = 1.0e+9
        sigma_T = 1.0e-16
        n = 1.0e+17
        ref = int(round(np.sqrt(2) / (32.0 * w * sigma_T**3 * n**2)))
        N = oc._calc_N_res(w, sigma_T, n)
        
        self.assertEqual(ref, N)
        
    def test__calc_n(self):
        box = ((0, 1), (2, 4), (4, 7))
        N = 300
        w = 100
        res = oc._calc_n(box, N, w)
        
        self.assertEqual(18.0, res)
        
    def test__is_inside(self):
        box = [[2.0, 3.0], [4.0, 6.0], [-1.0, 1.0]]
        box = np.array(box)
        position1 = np.array([2.5, 5.0, 0.0])
        position2 = np.array([0.0, 0.0, 0.0])
        position3 = np.array([2.5, 6.0, 0.5])
        
        self.assertTrue(oc._is_inside(position1, box))
        self.assertFalse(oc._is_inside(position2, box))
        self.assertTrue(oc._is_inside(position3, box))
        
    def test_sort(self):
        box = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
        positions1 = np.random.random((100, 3))
        positions2 = np.random.random((200, 3)) - np.ones((200, 3))*2
        positions = np.concatenate((positions2, positions1, positions2))
        permutations1 = np.array([i for i in range(len(positions2))])
        permutations2 = np.array([i for i in range(len(positions2), len(positions))])
        np.random.shuffle(permutations2)
        permutations = np.concatenate((permutations1, permutations2))
        offset = len(positions2)
        N = len(positions1) + len(positions2)
        count = 0
        
        for position in positions1:
            self.assertTrue(oc._is_inside(position, box))
            
        for position in positions2:
            self.assertFalse(oc._is_inside(position, box))
            
        for position in positions:
            if oc._is_inside(position, box):
                count += 1
                
        self.assertEqual(len(positions1), count)
        
        self.assertEqual(len(positions2), len(permutations1))
        self.assertEqual(len(positions1) + len(positions2), len(permutations2))
        self.assertEqual(len(positions), len(permutations))
        
        permutations, Nnew = oc._sort(permutations, box, positions, offset, N)
        
        self.assertEqual(Nnew, len(positions1))
        
        for i in range(offset, offset + Nnew):
            p = permutations[i]
            self.assertTrue(oc._is_inside(positions[p], box))
            
        for i in range(offset + Nnew, len(permutations)):
            p = permutations[i]
            self.assertFalse(oc._is_inside(positions[p], box))

    def test_sort2(self):
        box = np.array([[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]])
        N = 20
        Nh = 10
        Nm = 15
        Np = 5
        positions  = np.empty((N, 3))
        positions1 = np.empty((Np, 3))
        positions2 = np.empty((Np, 3))
        positions3 = np.empty((Nh, 3))
        permutations = np.array([i for i in range(N)])
        
        for i in range(Np):
            positions1[i] = np.random.random(3) - np.ones(3)
            
        for i in range(Np):
            positions2[i] = np.random.random(3)
            
        for i in range(Nh):
            positions3[i] = np.random.random(3) - np.ones(3)
            
        positions = np.concatenate((positions3, positions1, positions2))
        permutations, Nnew = oc._sort(permutations, box, positions, Nh, Nh)
        
        self.assertEqual(Nnew, Np)
        
        for i in range(Nh):
            p = permutations[i]
            pos = positions[p]
            a = not (pos[0] <= box[0][1])
            b = not (pos[0] >= box[0][0])
            
            c = not (pos[1] <= box[1][1])
            d = not (pos[1] >= box[1][0])
            
            e = not (pos[2] <= box[2][1])
            f = not (pos[2] >= box[2][0])
            
            self.assertTrue(a or b or c or d or e or f)
            
        for i in range(Nm, N):
            p = permutations[i]
            pos = positions[p]
            a = not (pos[0] <= box[0][1])
            b = not (pos[0] >= box[0][0])
            
            c = not (pos[1] <= box[1][1])
            d = not (pos[1] >= box[1][0])
            
            e = not (pos[2] <= box[2][1])
            f = not (pos[2] >= box[2][0])
            
            self.assertTrue(a or b or c or d or e or f)
            
        for i in range(Nh, Nm):
            p = permutations[i]
            pos = positions[p]
            a = (pos[0] <= box[0][1])
            b = (pos[0] >= box[0][0])
            
            c = (pos[1] <= box[1][1])
            d = (pos[1] >= box[1][0])
            
            e = (pos[2] <= box[2][1])
            f = (pos[2] >= box[2][0])
            
            self.assertTrue(a and b and c and d and e and f)


    def load_particles(self):
        positions = []
        
        with open("./test_data/particles.csv", "r", newline='') as csvfile:
            reader = csv.reader(csvfile, delimiter=',')
            for row in reader:
                positions.append([float(row[i]) for i in range(3)])
                
        return np.array(positions)

    def test_sort3(self):
        box = np.array([[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]])
        boxes = [] + oc._create_boxes(box)
        N = 100
        permutations = np.array([i for i in range(N)])            
        positions = self.load_particles()
            
        permutations , Np1 = oc._sort(permutations, boxes[0], positions, 0, N)
        permutations , Np2 = oc._sort(permutations, boxes[1], positions, Np1, N - Np1)
        permutations , Np3 = oc._sort(permutations, boxes[2], positions, Np1 + Np2, N - Np2 - Np1)
        #permutations , Np4 = oc._sort(permutations, boxes[3], positions, Np3, N - Np3)
        
        #permutations , Np5 = oc._sort(permutations, boxes[4], positions, Np4, N - Np4)
        #permutations , Np6 = oc._sort(permutations, boxes[5], positions, Np5, N - Np5)
        #permutations , Np7 = oc._sort(permutations, boxes[6], positions, Np6, N - Np6)
        #permutations , Np8 = oc._sort(permutations, boxes[7], positions, Np7, N - Np7)
        
        Nnew = np.zeros(8, dtype=int)
        
        for i in range(Np1):
        	p = permutations[i]
        	pos = positions[p]
        	if oc._is_inside(pos, boxes[0]):
        		Nnew[0] += 1
        		
        self.assertEqual(Nnew[0], Np1)
        
        for i in range(Np1, Np1 + Np2):
        	p = permutations[i]
        	pos = positions[p]
        	if oc._is_inside(pos, boxes[1]):
        		Nnew[1] += 1
        		
        self.assertEqual(Nnew[1], Np2)
        
        for i in range(Np1 + Np2, Np1 + Np2 + Np3):
        	p = permutations[i]
        	pos = positions[p]
        	if oc._is_inside(pos, boxes[2]):
        		Nnew[2] += 1
        		
        self.assertEqual(Nnew[2], Np3)

    def test__create_boxes(self):
        box_orig = np.array([[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]])
        boxes = oc._create_boxes(box_orig)
        
        self.assertEqual(8, len(boxes))
        V = 0.0
        
        for box in boxes:
            V += oc.get_V(box)
            
        self.assertEqual(oc.get_V(box_orig), V)
        
    def test__get_min_aspect_ratio(self):
        box = np.array([(0.0, 1.0), (0.0, 10.0), (0.0, 100.0)])
        axis1 = 0
        axis2 = 1
        axis3 = 2
        
        self.assertEqual(0.01, oc._get_min_aspect_ratio(box, axis1))
        self.assertEqual(0.1, oc._get_min_aspect_ratio(box, axis2))
        self.assertEqual(10.0, oc._get_min_aspect_ratio(box, axis3))
            
            
class TestOctreeOctree(unittest.TestCase):
    def test_build(self):
        positions = np.random.random((1000, 3))*2.0 - np.ones((1000, 3))
        octree = oc.Octree()
        octree.w = 1e+18
        
        octree.build(positions)
