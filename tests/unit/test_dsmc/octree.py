import unittest
import numpy as np
from dsmc import octree as oc
from dsmc import particles as part

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
        positions_subs = np.concatenate((positions1, positions2))
        np.random.shuffle(positions_subs)
        positions = np.concatenate((positions2, positions_subs))
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
        
        Nnew = oc._sort(box, positions, offset, N)
        
        self.assertEqual(Nnew, len(positions1))
        
        for i in range(offset, offset + Nnew):
            self.assertTrue(oc._is_inside(positions[i], box))
            
        for i in range(offset + Nnew, len(positions)):
            #print(offset, Nnew, i)
            self.assertFalse(oc._is_inside(positions[i], box))
