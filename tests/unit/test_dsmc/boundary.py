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