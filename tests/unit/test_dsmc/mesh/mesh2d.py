import numpy as np
import dsmc.mesh.mesh2d as msh
import unittest

class TestMesh2(unittest.TestCase):
    def test_constructor(self):
        mesh = msh.Mesh2d()
        
        self.assertEqual(1.0, mesh.cell_size)
        self.assertEqual(0.0, mesh.min1)
        self.assertEqual(0.0, mesh.min2)
        self.assertEqual(1, mesh.n_cells1)
        self.assertEqual(1, mesh.n_cells2)
        self.assertEqual(msh.Plane.XY, mesh.plane)