import unittest
from dsmc import particles as pa

class TestParticles(unittest.TestCase):
    def test_box_muller(self):
        T = 300
        result = pa.box_muller(T)
