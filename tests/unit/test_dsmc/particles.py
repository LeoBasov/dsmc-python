import unittest
from dsmc import particles as pa
import numpy as np

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

    def test_get_velocities(self):
        T = 300
        mass = 1.0e-26
        N = 1000
        u = np.zeros(3)
        velocities = pa.get_velocities(T, mass, N, u)

        self.assertEqual(N, len(velocities))

    def test_calc_temperature(self):
        T = 300
        mass = 1.0e-26
        N = 10000
        u = np.zeros(3)
        velocities = pa.get_velocities(T, mass, N, u)
        T_new = pa.calc_temperature(velocities, mass)
        diff = abs((T_new - T)/T)

        self.assertTrue(diff < 0.1)

    def test_calc_positions(self):
        x = (-1.0, 1.0)
        y = (2.0, 3.0)
        z = (-2.0, -1.0)
        X = np.array([x, y, z])
        N = 1000
        positions = pa.calc_positions(X, N)

        self.assertEqual(N, len(positions))

    def test_create_particles(self):
        x = (-1.0, 1.0)
        y = (2.0, 3.0)
        z = (-2.0, -1.0)
        X = np.array((x, y, z))
        N = 1000
        mass = 1.0e-26
        T = 300
        u = np.zeros(3)
        particles = pa.Particles()

        particles.create_particles(X, mass, T, N, u)

        self.assertEqual(N, len(particles.Pos))
        self.assertEqual(N, len(particles.Vel))
        self.assertEqual(N, particles.N)

        for i in range(N):
            self.assertTrue(particles.Pos[i][0] >= x[0] and particles.Pos[i][0] <= x[1])
            self.assertTrue(particles.Pos[i][1] >= y[0] and particles.Pos[i][1] <= y[1])
            self.assertTrue(particles.Pos[i][2] >= z[0] and particles.Pos[i][2] <= z[1])
            
        particles.create_particles(X, mass, T, N, u)

        self.assertEqual(2*N, len(particles.Pos))
        self.assertEqual(2*N, len(particles.Vel))
        self.assertEqual(2*N, particles.N)
        
        for i in range(2*N):
            self.assertTrue(particles.Pos[i][0] >= x[0] and particles.Pos[i][0] <= x[1])
            self.assertTrue(particles.Pos[i][1] >= y[0] and particles.Pos[i][1] <= y[1])
            self.assertTrue(particles.Pos[i][2] >= z[0] and particles.Pos[i][2] <= z[1])
            
    def test_calc_vp(self):
        T = 300
        mass = 1e-26
        vp = pa.calc_vp(T, mass)
        
        self.assertEqual(np.sqrt(2.0*pa.kb*T/mass), vp)
