import unittest
from dsmc import particles as pa

class TestParticles(unittest.TestCase):
    def test_box_muller(self):
        T = 300
        result = pa.box_muller(T)

    def test_x2velocity(self):
        x = 1.0
        mass = 1.0e-26
        result = pa.x2velocity(x, mass)

    def test_get_vel(self):
        T = 300
        mass = 1.0e-26
        result = pa.get_vel(T, mass)
