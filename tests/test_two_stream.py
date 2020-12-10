from collections import namedtuple
from random import uniform
from unittest import main, TestCase

from numpy import exp, seterr, zeros

from pyrad.solvers.two_stream import adding, delta_eddington


Optics = namedtuple("Optics", ("tau", "omega", "g"))


class TestTwoStream(TestCase):

    def test_simple(self):
        num_layers = 5
        R_direct, T_direct, T_pure = zeros(num_layers), zeros(num_layers), zeros(num_layers)
        R_diffuse, T_diffuse = zeros(num_layers), zeros(num_layers)
        zenith = [uniform(0.25, 0.75), 0.5]
        for i in range(num_layers):
            optics = Optics(tau=uniform(0.1, 1.), omega=uniform(0.5, 0.9),
                            g=uniform(0.2, 0.9))
            R_direct[i], T_direct[i], T_pure[i] = delta_eddington(optics.tau, optics.omega,
                                                                  optics.g, zenith[0])
            R_diffuse[i], T_diffuse[i], _ = delta_eddington(optics.tau, optics.omega,
                                                            optics.g, zenith[-1])
        R, T = adding(R_direct, R_diffuse, T_direct, T_diffuse, T_pure, 0.6, 0.6)

    def test_no_gas(self):
        optics = Optics(tau=0., omega=0.85, g=0.85)
        zenith = 0.75
        R, T, T_pure = delta_eddington(optics.tau, optics.omega, optics.g, zenith)
        self.assertEqual(R, 0.)
        self.assertEqual(T, 1.)
        self.assertEqual(T_pure, T)

    def test_no_scattering(self):
        optics = Optics(tau=1., omega=0., g=0.85)
        zenith = 0.75
        R, T, T_pure = delta_eddington(optics.tau, optics.omega, optics.g, zenith)
        self.assertEqual(R, 0.)
        self.assertAlmostEqual(T, exp(-1.*optics.tau/zenith))
        self.assertEqual(T_pure, T)

    def test_conservative_scattering(self):
        optics = Optics(tau=1., omega=1., g=0.85)
        zenith = 0.75
        R, T, T_pure = delta_eddington(optics.tau, optics.omega, optics.g, zenith)
        self.assertAlmostEqual(R + T, 1.)


if __name__ == "__main__":
    seterr(divide="raise", over="raise", invalid="raise")
    main()
