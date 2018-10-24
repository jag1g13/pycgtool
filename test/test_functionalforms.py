import unittest

import numpy as np

from pycgtool.functionalforms import FunctionalForms, FunctionalForm
from pycgtool.util import circular_mean, circular_variance


class FunctionalFormTest(unittest.TestCase):
    def test_functional_form(self):
        funcs = FunctionalForms()
        harmonic_form = funcs.Harmonic
        self.assertIsNotNone(harmonic_form)
        cos_harmonic_form = funcs.CosHarmonic
        self.assertIsNotNone(cos_harmonic_form)

    def test_functional_form_new(self):
        class TestFunc(FunctionalForm):
            gromacs_type_ids = (None, None, None)

            @staticmethod
            def eqm(values, temp):
                return "TestResultEqm"

            @staticmethod
            def fconst(values, temp):
                return "TestResultFconst"

        funcs = FunctionalForms()
        self.assertIn("TestFunc", funcs)
        self.assertEqual("TestResultEqm", funcs.TestFunc.eqm(None, None))
        self.assertEqual("TestResultFconst", funcs.TestFunc.fconst(None, None))

    def test_inject_mean_function(self):
        """
        Test injecting mean function into Boltzmann Inversion
        """
        funcs = FunctionalForms()
        vals = [0.25 * np.pi, 1.75 * np.pi]

        func = funcs["Harmonic"]()

        np.testing.assert_allclose(
            func.eqm(vals, None), np.pi
        )

        func = funcs["Harmonic"](
            mean_func=circular_mean
        )

        np.testing.assert_allclose(
            func.eqm(vals, None), 0,
            atol=1e-6
        )

    def test_inject_variance_function(self):
        """
        Test injecting variance function into Boltzmann Inversion
        """
        funcs = FunctionalForms()

        func = funcs["Harmonic"](
            variance_func=circular_variance
        )

        vals = [0, 0.5 * np.pi]
        ref = func.fconst(vals, 300)

        vals = [0.25 * np.pi, 1.75 * np.pi]
        np.testing.assert_allclose(
            func.fconst(vals, 300), ref,
            atol=0
        )


if __name__ == '__main__':
    unittest.main()
