import unittest
from dsmc import dsmc

class TestDSMC(unittest.TestCase):
    def test_test(self):
        self.assertEqual(1.0, dsmc.test())
