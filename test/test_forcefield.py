import unittest
import os
import collections

from pycgtool.forcefield import ForceField

DummyBond = collections.namedtuple("DummyBond", ["atoms", "eqm", "fconst"])


class ForceFieldTest(unittest.TestCase):
    def test_create(self):
        name = "test"
        dirname = "fftest.ff"

        ForceField(name)
        self.assertTrue(os.path.exists(dirname))
        self.assertTrue(os.path.isdir(dirname))
        ForceField(name)

    def test_bond_section(self):
        bonds = [
            DummyBond(["a", "b"], 1, 100),
            DummyBond(["b", "c"], 2, 50)
        ]

        expected = ["  [ dummy ]",
                    "       a    b      1.00000    100.00000",
                    "       b    c      2.00000     50.00000"
        ]

        self.assertListEqual(expected, ForceField.bond_section(bonds, "dummy"))

    def test_r2b(self):
        nter = {"a", "b"}
        cter = {"a", "c"}
        bter = {"a", "d"}

        expected = [
            "; rtp residue to rtp building block table",
            ";     main  N-ter C-ter 2-ter",
            "a     a     Na    Ca    2a   ",
            "b     b     Nb    -     -    ",
            "c     c     -     Cc    -    ",
            "d     d     -     -     2d   "
        ]

        output = ForceField.write_r2b(nter, cter, bter)
        self.assertListEqual(expected[:2], output[:2])
        self.assertListEqual(expected[2:], output[2:])


if __name__ == '__main__':
    unittest.main()
