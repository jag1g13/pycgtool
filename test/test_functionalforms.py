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
            @staticmethod
            def __call__(mean, var, temp):
                return "TestResult"

        funcs = FunctionalForms()
        self.assertIn("TestFunc", funcs)
        self.assertEqual("TestResult", funcs.TestFunc(None, None, None))

if __name__ == '__main__':
    unittest.main()
