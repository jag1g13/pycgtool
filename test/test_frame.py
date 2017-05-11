import unittest
import filecmp
import os
import logging

import numpy as np

from pycgtool.frame import Atom, Residue, Frame
from pycgtool.frame import FrameReaderSimpleTraj, FrameReaderMDAnalysis, FrameReader

try:
    import mdtraj
    mdtraj_present = True
except ImportError:
    mdtraj_present = False

try:
    import MDAnalysis
    mdanalysis_present = True
except ImportError:
    mdanalysis_present = False


class AtomTest(unittest.TestCase):
    def test_atom_create(self):
        atom = Atom(name="Name", num=0, type="Type")
        self.assertEqual("Name", atom.name)
        self.assertEqual(0, atom.num)
        self.assertEqual("Type", atom.type)

    def test_atom_add_missing_data(self):
        atom1 = Atom("Name1", 0, type="Type")
        atom2 = Atom("Name2", 0, mass=1)

        with self.assertRaises(AssertionError):
            atom1.add_missing_data(atom2)

        atom2 = Atom("Name1", 0, mass=1)
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

    def helper_read_gro(self, frame):
        self.assertEqual(221, len(frame.residues))
        self.assertEqual("SOL", frame.residues[0].name)
        self.assertEqual(3, len(frame.residues[0].atoms))
        self.assertEqual("OW", frame.residues[0].atoms[0].name)
        np.testing.assert_allclose(np.array([0.696, 1.33, 1.211]),
                                   frame.residues[0].atoms[0].coords)

    def test_frame_simpletraj_read_gro(self):
        frame = Frame("test/data/water.gro", xtc_reader="simpletraj")

        self.helper_read_gro(frame)

    # MDTRAJ changes water name to HOH
    @unittest.expectedFailure
    @unittest.skipIf(not mdtraj_present, "MDTRAJ or Scipy not present")
    def test_frame_mdtraj_read_gro(self):
        logging.disable(logging.WARNING)
        frame = Frame("test/data/water.gro", xtc_reader="mdtraj")
        logging.disable(logging.NOTSET)

        self.helper_read_gro(frame)

    @unittest.skipIf(not mdanalysis_present, "MDAnalysis not present")
    def test_frame_mdanalysis_read_gro(self):
        reader = FrameReaderMDAnalysis("test/data/water.gro")
        frame = Frame.instance_from_reader(reader)

        self.helper_read_gro(frame)

    def test_frame_output_gro(self):
        frame = Frame("test/data/water.gro")
        frame.output("water-out.gro", format="gro")
        self.assertTrue(filecmp.cmp("test/data/water.gro", "water-out.gro"))
        os.remove("water-out.gro")

    def test_frame_read_xtc_simpletraj_numframes(self):
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc",
                      xtc_reader="simpletraj")
        self.assertEqual(11, frame.numframes)

    @unittest.skipIf(not mdtraj_present, "MDTRAJ or Scipy not present")
    def test_frame_read_xtc_mdtraj_numframes(self):
        logging.disable(logging.WARNING)
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc",
                      xtc_reader="mdtraj")
        logging.disable(logging.NOTSET)
        self.assertEqual(11, frame.numframes)

    def helper_read_xtc(self, frame):
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

    def test_frame_simpletraj_read_xtc(self):
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc",
                      xtc_reader="simpletraj")
        self.helper_read_xtc(frame)

    @unittest.skipIf(not mdtraj_present, "MDTRAJ or Scipy not present")
    def test_frame_mdtraj_read_xtc(self):
        logging.disable(logging.WARNING)
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc",
                      xtc_reader="mdtraj")
        logging.disable(logging.NOTSET)

        self.helper_read_xtc(frame)

    @unittest.skipIf(not mdanalysis_present, "MDAnalysis not present")
    def test_frame_mdanalysis_read_xtc(self):
        reader = FrameReaderMDAnalysis("test/data/water.gro", "test/data/water.xtc")
        frame = Frame.instance_from_reader(reader)

        self.helper_read_xtc(frame)

    @unittest.skipIf(not mdtraj_present, "MDTRAJ or Scipy not present")
    def test_frame_write_xtc_mdtraj(self):
        try:
            os.remove("water_test2.xtc")
        except IOError:
            pass

        logging.disable(logging.WARNING)
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc",
                      xtc_reader="mdtraj")
        logging.disable(logging.NOTSET)

        while frame.next_frame():
            frame.write_xtc("water_test2.xtc")

    def test_frame_instance_from_reader(self):
        reader = FrameReaderSimpleTraj("test/data/water.gro")
        frame = Frame.instance_from_reader(reader)

        self.helper_read_gro(frame)

    def test_frame_instance_from_reader_dummy(self):
        class DummyReader(FrameReader):
            def _initialise_frame(self, frame):
                frame.dummy_reader = True

            def _read_frame_number(self, number):
                return number * 10, [], None

        reader = DummyReader(None)
        frame = Frame.instance_from_reader(reader)
        self.assertTrue(frame.dummy_reader)

        frame.next_frame()
        self.assertEqual(frame.number, 0)
        self.assertEqual(frame.time, 0)
        self.assertIsNone(frame.box)

        frame.next_frame()
        self.assertEqual(frame.number, 1)
        self.assertEqual(frame.time, 10)


if __name__ == '__main__':
    unittest.main()
