import unittest
import os

from pycgtool.forcefield import ForceField


class ForceFieldTest(unittest.TestCase):
    def test_create(self):
        name = "test"
        dirname = "fftest.ff"
        ForceField(name)
        self.assertTrue(os.path.exists(dirname))
        self.assertTrue(os.path.isdir(dirname))
        ForceField(name)

if __name__ == '__main__':
    unittest.main()
