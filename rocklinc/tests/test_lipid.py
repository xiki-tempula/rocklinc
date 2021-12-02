from rocklinc.lipid import StepLipid, CurveLipid
from rocklinc.lipid import plot_dx
import numpy as np
import matplotlib.pyplot as plt
import pytest

class TestStepLipid():
    @staticmethod
    @pytest.fixture(scope='class')
    def lipid():
        lipid = StepLipid([80, 80, 80], [49, 49, 49])
        lipid.write_dx_input([6.22602388e-04, 18.5969598, 1.98288109,
                                         24.1998191, 1.91347013],
                              [13.80515, 2, 80])

        lipid.write_apbs_input()
        lipid.run_APBS(apbs_exe='apbs')
        charge, diel, esp = lipid.read_APBS()
        fig = plot_dx(np.linspace(0, lipid.dim[-1].magnitude, lipid.grid[-1]), charge, diel, esp)
        return lipid, fig

    def test_fig(self, lipid):
        lipid, fig = lipid
        assert isinstance(fig, plt.figure)

class TestCurveLipid():
    @staticmethod
    @pytest.fixture(scope='class')
    def lipid():
        lipid = CurveLipid([80, 80, 80], [49, 49, 49])
        lipid.write_dx_input([24.1998191, 1.28409182, 3.29937546e-03,
                              18.5969598, 2.41286946, 5.70594706e-03,
                              16.0770291, 2.47034776, 4.17737240e-03,
                              6.59972225, 0.90515386, 4.25321040e-04,
                              2.95472778, 1.17224618e-04],
                              [0.36626031, 18.87278382, 76.49060424, 4.83151723])

        lipid.write_apbs_input()
        lipid.run_APBS(apbs_exe='apbs')
        charge, diel, esp = lipid.read_APBS()
        fig = plot_dx(np.linspace(0, lipid.dim[-1].magnitude, lipid.grid[-1]), charge, diel, esp)
        return lipid, fig
    
    def test_fig(self, lipid):
        lipid, fig = lipid
        assert isinstance(fig, plt.figure)
