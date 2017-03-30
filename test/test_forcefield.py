import unittest
import os
import collections

from pycgtool.forcefield import ForceField

DummyBMap = collections.namedtuple("DummyBMap", ["name", "type", "charge"])


class DummyBond:
    def __init__(self, atoms, eqm, fconst):
        self.atoms = atoms
        self.eqm = eqm
        self.fconst = fconst

    def __iter__(self):
        return iter(self.atoms)


class DummyBondSet:
    def __init__(self, bonds, name):
        self.bonds = bonds
        self.name = name

    def __contains__(self, item):
        return self.name == item

    def bonds_with_natoms(self, n):
        return [bond for bond in self.bonds if len(bond.atoms) == n]

    def get_bond_lengths(self, *args, **kwargs):
        return self.bonds_with_natoms(2)

    def get_bond_angles(self, *args, **kwargs):
        return self.bonds_with_natoms(3)

    def get_bond_dihedrals(self, *args, **kwargs):
        return self.bonds_with_natoms(4)


class ForceFieldTest(unittest.TestCase):
    def setUp(self):
        self.bonds = [
            DummyBond(["a",  "b"], 1, 100),
            DummyBond(["b",  "c"], 2,  50),
            DummyBond(["c", "+a"], 3,  20),
        ]

        self.mapping = {"Dummy": [
                DummyBMap("a", "a1",  1),
                DummyBMap("b", "b2", -1),
                DummyBMap("c", "c3",  0),
        ]}

        self.bondset = DummyBondSet(self.bonds, "Dummy")

    def test_create(self):
        name = "test"
        dirname = "fftest.ff"

        ForceField(name)
        self.assertTrue(os.path.exists(dirname))
        self.assertTrue(os.path.isdir(dirname))
        ForceField(name)

    def test_bond_section(self):
        expected = [
            "  [ section ]",
            "       a    b      1.00000    100.00000",
            "       b    c      2.00000     50.00000",
            "       c   +a      3.00000     20.00000"
        ]

        self.assertListEqual(expected, ForceField.bond_section(self.bonds, "section"))

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
        self.assertListEqual(expected, output)

    def test_rtp(self):
        expected = [
            "[ bondedtypes ]",
            "   1   1   1   1   1   1   0   0",
            "[ Dummy ]",
            "  [ atoms ]",
            "       a   a1 1.000000    0",
            "       b   b2 -1.000000    0",
            "       c   c3 0.000000    0",
            "  [ bonds ]",
            "       a    b      1.00000    100.00000",
            "       b    c      2.00000     50.00000",
            "       c   +a      3.00000     20.00000",
            "[ CDummy ]",
            "  [ atoms ]",
            "       a   a1 1.000000    0",
            "       b   b2 -1.000000    0",
            "       c   c3 0.000000    0",
            "  [ bonds ]",
            "       a    b      1.00000    100.00000",
            "       b    c      2.00000     50.00000"
        ]

        output, nter, cter, bter = ForceField.write_rtp(self.mapping, self.bondset)

        self.maxDiff = None
        self.assertListEqual(expected, output)

        self.assertFalse(nter)
        self.assertEqual({"Dummy"}, cter)
        self.assertFalse(bter)


if __name__ == '__main__':
    unittest.main()
