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
        
    def test_build(self):
        particles = part.Particles()
        tree = oc.Octree()
        
        X = ((-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0))
        mass = 1.0e-26
        T = 300.0
        N = 10000
        
        particles.create_particles(X, mass, T, N)
        tree.build(particles.Pos)
