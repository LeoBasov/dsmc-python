import numpy as np
import dsmc.boundary as bo
import unittest

class TestCommon(unittest.TestCase):
    def test__check_if_parallel(self):
        v1 = np.array([1.0, 0.0, 0.0])
        v2 = np.array([1.0, 0.0, 0.0])
        v3 = np.array([0.0, 1.0, 1.0])
        v4 = np.array([1.0e-6, 10.0, 10.0])
        
        self.assertTrue(bo._check_if_parallel(v1, v2))
        self.assertFalse(bo._check_if_parallel(v1, v3))
        self.assertFalse(bo._check_if_parallel(v1, v4))
        
    def test__intersect(self):
        l0 = np.array([-1.0, -1.0, 0.0])
        l1 = np.array([1.0, 1.0, 0.0])
        
        p0 = np.array([0.0, -1.0, -1.0])
        p1 = np.array([0.0, -1.0, 1.0])
        p2 = np.array([0.0, 1.0, 1.0])
        
        intersected, n_l, n_r, t = bo._intersect(l0, l1, p0, p1, p2)
        
        self.assertTrue(intersected)
        
        self.assertEqual(2.0, n_l[0])
        self.assertEqual(2.0, n_l[1])
        self.assertEqual(0.0, n_l[2])
        
        self.assertEqual(-4.0, n_r[0])
        self.assertEqual(0.0, n_r[1])
        self.assertEqual(0.0, n_r[2])
        
        self.assertEqual(0.5, t)