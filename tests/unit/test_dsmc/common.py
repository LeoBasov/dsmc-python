import numpy as np
import dsmc.common as com
import unittest

class TestCommon(unittest.TestCase):
    def test_pass(self):
        box = np.array([[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]])
        V= com.get_V(box)