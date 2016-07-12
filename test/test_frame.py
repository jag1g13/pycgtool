import unittest
import filecmp
import os

import numpy as np

from pycgtool.frame import Atom, Residue, Frame

try:
    import mdtraj
    mdtraj_present = True
except ImportError:
    mdtraj_present = False


class AtomTest(unittest.TestCase):
    def test_atom_create(self):
        atom = Atom(name="Name", num=0, type="Type")
        self.assertEqual("Name", atom.name)
        self.assertEqual(0, atom.num)
        self.assertEqual("Type", atom.type)

    def test_atom_add_missing_data(self):
        atom1 = Atom(name="Name", num=0, type="Type")
        atom2 = Atom(mass=1)

        with self.assertRaises(AssertionError):
            atom1.add_missing_data(atom2)

        atom2 = Atom(name="Name", num=0, mass=1)
        atom1.add_missing_data(atom2)
        self.assertEqual(1, atom1.mass)


class ResidueTest(unittest.TestCase):
    def test_residue_create(self):
        residue = Residue(name="Resname")
        self.assertEqual("Resname", residue.name)

    def test_residue_add_atoms(self):
        atom = Atom(name="Name", num=0, type="Type")
        residue = Residue()
        residue.add_atom(atom)
        self.assertEqual(atom, residue.atoms[0])
        self.assertTrue(atom is residue.atoms[0])


class FrameTest(unittest.TestCase):
    def test_frame_create(self):
        Frame()

    def test_frame_add_residue(self):
        residue = Residue()
        frame = Frame()
        frame.add_residue(residue)
        self.assertEqual(residue, frame.residues[0])
        self.assertTrue(residue is frame.residues[0])

    def test_frame_read_gro(self):
        frame = Frame("test/data/water.gro")
        self.assertEqual(221, len(frame.residues))
        self.assertEqual("SOL", frame.residues[0].name)
        self.assertEqual(3, len(frame.residues[0].atoms))
        self.assertEqual("OW", frame.residues[0].atoms[0].name)
        np.testing.assert_allclose(np.array([0.696, 1.33, 1.211]),
                                   frame.residues[0].atoms[0].coords)

    def test_frame_output_gro(self):
        frame = Frame("test/data/water.gro")
        frame.output("water-out.gro", format="gro")
        self.assertTrue(filecmp.cmp("test/data/water.gro", "water-out.gro"))
        os.remove("water-out.gro")

    def test_frame_read_xtc_simpletraj_numframes(self):
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc",
                      xtc_reader="simpletraj")
        self.assertEqual(12, frame.numframes)

    def test_frame_read_xtc_simpletraj(self):
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc",
                      xtc_reader="simpletraj")
        self.assertEqual(663, frame.natoms)
        # These are the coordinates from the gro file
        np.testing.assert_allclose(np.array([0.696, 1.33, 1.211]),
                                   frame.residues[0].atoms[0].coords)
        np.testing.assert_allclose(np.array([1.89868, 1.89868, 1.89868]),
                                   frame.box)

        frame.next_frame()
        # These coordinates are from the xtc file
        np.testing.assert_allclose(np.array([1.176, 1.152, 1.586]),
                                   frame.residues[0].atoms[0].coords)
        np.testing.assert_allclose(np.array([1.9052, 1.9052, 1.9052]),
                                   frame.box)
        frame.next_frame()
        np.testing.assert_allclose(np.array([1.122, 1.130, 1.534]),
                                   frame.residues[0].atoms[0].coords)
        np.testing.assert_allclose(np.array([1.90325272, 1.90325272, 1.90325272]),
                                   frame.box)

    @unittest.skipIf(not mdtraj_present, "MDTRAJ or Scipy not present")
    def test_frame_read_xtc_mdtraj_numframes(self):
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc",
                      xtc_reader="mdtraj")
        self.assertEqual(12, frame.numframes)

    @unittest.skipIf(not mdtraj_present, "MDTRAJ or Scipy not present")
    def test_frame_read_xtc_mdtraj(self):
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc",
                      xtc_reader="mdtraj")
        self.assertEqual(663, frame.natoms)
        # These are the coordinates from the gro file
        np.testing.assert_allclose(np.array([0.696, 1.33, 1.211]),
                                   frame.residues[0].atoms[0].coords)
        np.testing.assert_allclose(np.array([1.89868, 1.89868, 1.89868]),
                                   frame.box)
        frame.next_frame()
        # These coordinates are from the xtc file
        np.testing.assert_allclose(np.array([1.176, 1.152, 1.586]),
                                   frame.residues[0].atoms[0].coords)
        np.testing.assert_allclose(np.array([1.9052, 1.9052, 1.9052]),
                                   frame.box)
        frame.next_frame()
        np.testing.assert_allclose(np.array([1.122, 1.130, 1.534]),
                                   frame.residues[0].atoms[0].coords)
        np.testing.assert_allclose(np.array([1.90325272, 1.90325272, 1.90325272]),
                                   frame.box)

    @unittest.skipIf(not mdtraj_present, "MDTRAJ or Scipy not present")
    def test_frame_write_xtc_mdtraj(self):
        try:
            os.remove("water_test2.xtc")
        except IOError:
            pass
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc",
                      xtc_reader="mdtraj")
        while frame.next_frame():
            frame.write_xtc("water_test2.xtc")


if __name__ == '__main__':
    unittest.main()
