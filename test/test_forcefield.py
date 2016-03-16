import unittest
import os

from pycgtool.forcefield import ForceField


class ForceFieldTest(unittest.TestCase):
    def test_create(self):
        name = "fftest.ff"
        ForceField(name)
        self.assertTrue(os.path.exists(name))
        self.assertTrue(os.path.isdir(name))
        ForceField(name)

if __name__ == '__main__':
    unittest.main()
