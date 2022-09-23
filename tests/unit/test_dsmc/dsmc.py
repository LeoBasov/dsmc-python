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
        
    def test__calc_post_col_vels(self):
        velocity1 : np.ndarray = np.array([1.0, 2.0, 3.0])
        velocity2 : np.ndarray = np.array([1.0, 2.0, 3.0])
        mass1 : float = 1.0
        mass2 : float = 1.0
        rel_vel_module : float = 1.0
        rand_number1 : float = 1.0
        rand_number2 : float = 1.0
        
        res = ds._calc_post_col_vels(velocity1, velocity2, mass1, mass2, rel_vel_module, rand_number1, rand_number2)
        
        self.assertEqual(1.5, res[0][0])
        self.assertEqual(2.0, res[0][1])
        self.assertEqual(3.0, res[0][2])
        
        self.assertEqual(0.5, res[1][0])
        self.assertEqual(2.0, res[1][1])
        self.assertEqual(3.0, res[1][2])
