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