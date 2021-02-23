from unittest import main, TestCase

from numpy import arange

from pyrad.lbl.hitran.collision_induced_absorption import HitranCIA


class TestCollisionInducedAbsorption(TestCase):

    def test_collision_induced_absorption(self):
        interactions = [["CO2", "CO2"], ["N2", "N2"], ["O2", "N2"], ["O2", "O2"]]
        grid = arange(1., 3250., 0.01)
        for interaction in interactions:
            molecule = HitranCIA(*interaction)
            for temperature in [211.11, 235.45, 250.23, 270.45, 290.11, 300.1]:
                molecule.absorption_coefficient(temperature, grid)


if __name__ == "__main__":
    main()
