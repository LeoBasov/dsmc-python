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

    def test_get_velocities(self):
        T = 300
        mass = 1.0e-26
        N = 1000
        velocities = pa.get_velocities(T, mass, N)

        self.assertEqual(N, len(velocities))

    def test_calc_temperature(self):
        T = 300
        mass = 1.0e-26
        N = 10000
        velocities = pa.get_velocities(T, mass, N)
        T_new = pa.calc_temperature(velocities, mass)
        diff = abs((T_new - T)/T)

        self.assertTrue(diff < 0.1)

    def test_calc_positions(self):
        x = (-1.0, 1.0)
        y = (2.0, 3.0)
        z = (-2.0, -1.0)
        N = 1000
        positions = pa.calc_positions(x, y, z, N)

        self.assertEqual(N, len(positions))

    def test_create_particles(self):
        x = (-1.0, 1.0)
        y = (2.0, 3.0)
        z = (-2.0, -1.0)
        X = (x, y, z)
        N = 1000
        mass = 1.0e-26
        T = 300
        particles = pa.Particles()

        particles.create_particles(X, mass, T, N)

        self.assertEqual(N, len(particles.positions))
        self.assertEqual(N, len(particles.velocities))
        self.assertEqual(N, particles.N)
