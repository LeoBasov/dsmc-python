import dsmc.mesh.mesh2d as msh
import unittest
import numpy as np

class TestMesh2(unittest.TestCase):
    def test_constructor(self):
        mesh = msh.Mesh2d()
        
        self.assertEqual(1.0, mesh.cell_size)
        self.assertEqual(0.0, mesh.min1)
        self.assertEqual(0.0, mesh.min2)
        self.assertEqual(1, mesh.n_cells1)
        self.assertEqual(1, mesh.n_cells2)
        self.assertEqual(msh.Plane.XY, mesh.plane)
        
    def test__get_cell_id1(self):
        val1_f = -1.0
        val1_t = 0.1
        val2_f = -1.0
        val2_t = 0.1
        n_cells1 = 10
        n_cells2 = 10
        min1 = -0.5
        min2 = -0.5
        cell_size = 0.1
        
        res1_f = msh._get_cell_id(val1_t, val2_f, n_cells1, n_cells2, min1, min2, cell_size)
        res2_f = msh._get_cell_id(val1_f, val2_t, n_cells1, n_cells2, min1, min2, cell_size)
        res3_t = msh._get_cell_id(val1_t, val2_t, n_cells1, n_cells2, min1, min2, cell_size)
        
        self.assertFalse(res1_f[0])
        self.assertFalse(res2_f[0])
        self.assertTrue(res3_t[0])
        
    def test__get_cell_id2(self):
        val1 = (-0.45, -0.5)
        val2 = (-0.35, -0.35)
        n_cells1 = 10
        n_cells2 = 10
        min1 = -0.5
        min2 = -0.5
        cell_size = 0.1
        
        res1 = msh._get_cell_id(val1[0], val1[1], n_cells1, n_cells2, min1, min2, cell_size)
        res2 = msh._get_cell_id(val2[0], val2[1], n_cells1, n_cells2, min1, min2, cell_size)
        
        self.assertTrue(res1[0])
        self.assertTrue(res2[0])
        
        self.assertEqual(0, res1[1])
        self.assertEqual(11, res2[1])
        
    def test_get_cell_id(self):
        mesh = msh.Mesh2d()
        
        mesh.n_cells1 = 10
        mesh.n_cells2 = 10
        mesh.min1 = -0.5
        mesh.min2 = -0.5
        mesh.cell_size = 0.1
        
        pos = np.array([-0.45, -0.35, -0.25])
        
        # XY
        mesh.plane = msh.Plane.XY
        res = mesh.get_cell_id(pos)
        
        self.assertTrue(res[0])
        self.assertEqual(10, res[1])
        
        # YZ
        mesh.plane = msh.Plane.YZ
        res = mesh.get_cell_id(pos)
        
        self.assertTrue(res[0])
        self.assertEqual(21, res[1])
        
        # XZ
        mesh.plane = msh.Plane.XZ
        res = mesh.get_cell_id(pos)
        
        self.assertTrue(res[0])
        self.assertEqual(20, res[1])