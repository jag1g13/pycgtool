import unittest
import filecmp
import os

import numpy as np

from pycgtool.frame import Atom, Residue, Frame


class AtomTest(unittest.TestCase):
    def test_atom_create(self):
        atom = Atom(name="Name", num=0, type="Type")
        self.assertEqual("Name", atom.name)
        self.assertEqual(0, atom.num)
        self.assertEqual("Type", atom.type)


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
        frame = Frame()

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

    def test_frame_read_xtc(self):
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc")
        # These are the coordinates from the gro file
        np.testing.assert_allclose(np.array([0.696, 1.33, 1.211]),
                                   frame.residues[0].atoms[0].coords)
        frame.next_frame()
        # These coordinates are from the xtc file
        np.testing.assert_allclose(np.array([1.176, 1.152, 1.586]),
                                   frame.residues[0].atoms[0].coords)
        frame.next_frame()
        np.testing.assert_allclose(np.array([1.122, 1.130, 1.534]),
                                   frame.residues[0].atoms[0].coords)


if __name__ == '__main__':
    unittest.main()
