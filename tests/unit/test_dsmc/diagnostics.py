import numpy as np
import dsmc.diagnostics as dia
import unittest

class TestDiagnostics(unittest.TestCase):
    def test_sort_bin(self):
        positions = np.array(([1, 2, 3], [9, 8, 7], [10, 11, 12], [4, 5, 6]))
        Nbins = 4
        
        bins1, box = dia.sort_bin(positions, 0, Nbins)
        
        self.assertEqual(Nbins, len(bins1))
        
        for b in bins1:
            self.assertEqual(1, len(b))
            
        self.assertEqual(0, bins1[0][0])
        self.assertEqual(3, bins1[1][0])
        self.assertEqual(1, bins1[2][0])
        self.assertEqual(2, bins1[3][0])
