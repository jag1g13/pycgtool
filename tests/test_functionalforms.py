import unittest

import numpy as np

from pycgtool.functionalforms import FunctionalForm, get_functional_forms
from pycgtool.util import circular_mean, circular_variance


class FunctionalFormTest(unittest.TestCase):
    def test_functional_form(self):
        funcs = get_functional_forms()

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

        funcs = get_functional_forms()
        test_func = funcs["TestFunc"].value

        self.assertEqual("TestResultEqm", test_func.eqm(None, None))
        self.assertEqual("TestResultFconst", test_func.fconst(None, None))

    def test_inject_mean_function(self):
        """
        Test injecting mean function into Boltzmann Inversion
        """
        funcs = get_functional_forms()
        vals = [0.25 * np.pi, 1.75 * np.pi]

        func = funcs["Harmonic"].value()

        np.testing.assert_allclose(func.eqm(vals, None), np.pi)

        func = funcs["Harmonic"].value(mean_func=circular_mean)

        np.testing.assert_allclose(func.eqm(vals, None), 0, atol=1e-6)

    def test_inject_variance_function(self):
        """
        Test injecting variance function into Boltzmann Inversion
        """
        funcs = get_functional_forms()

        func = funcs["Harmonic"].value(variance_func=circular_variance)

        vals = [0, 0.5 * np.pi]
        ref = func.fconst(vals, 300)

        vals = [0.25 * np.pi, 1.75 * np.pi]
        np.testing.assert_allclose(func.fconst(vals, 300), ref, atol=0)


if __name__ == "__main__":
    unittest.main()
