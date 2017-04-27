import unittest

from pycgtool.functionalforms import FunctionalForms, FunctionalForm


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

if __name__ == '__main__':
    unittest.main()
