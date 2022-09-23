import unittest
import numpy as np
from dsmc import dsmc as ds

class TestDSMC(unittest.TestCase):
    def test_Constructor(self):
        ds.DSMC()
        
    def test__calc_prob(self):
        vel1 : np.ndarray = np.array([1.0, 2.0, 3.0])
        vel2 : np.ndarray = np.array([1.0, 2.0, 3.0])
        sigma_T : float = 1.0
        Vc : float = 1.0
        dt : float = 1.0
        w : float = 1.0
        N : int = 1
        
        res = ds._calc_prob(vel1, vel2, sigma_T, Vc, dt, w, N)
        
        self.assertEqual(np.linalg.norm(vel1 - vel2), res)
