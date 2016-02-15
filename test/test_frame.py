import unittest
import filecmp

import numpy as np

from pycgtool.frame import *


class MappingTest(unittest.TestCase):
    def test_mapping_create(self):
        mapping = Mapping("test/data/water.map")
        self.assertEqual(1, len(mapping))
        self.assertTrue("SOL" in mapping)
        self.assertEqual(1, len(mapping["SOL"]))
        self.assertEqual(3, len(mapping["SOL"][0].atoms))
        self.assertEqual("OW", mapping["SOL"][0].atoms[0])

    def test_mapping_apply(self):
        mapping = Mapping("test/data/water.map")
        frame = Frame("test/data/water.gro")
        cgframe = mapping.apply(frame)
        self.assertEqual(len(frame), len(cgframe))
        cgframe.output_gro("water-cg.gro")
        self.assertTrue(filecmp.cmp("test/data/water-cg.gro", "water-cg.gro"))


class MeasureTest(unittest.TestCase):
    def test_measure_create(self):
        measure = Measure("test/data/sugar.bnds")
        self.assertEqual(2, len(measure))
        self.assertTrue("SOL" in measure)
        self.assertTrue("ALLA" in measure)
        self.assertEqual(0, len(measure["SOL"]))
        self.assertEqual(18, len(measure["ALLA"]))

    @unittest.expectedFailure
    def test_measure_apply(self):
        measure = Measure("test/data/sugar.bnds")
        frame = Frame("test/data/sugar-cg.gro")
        measure.apply(frame)
        # First six are bond lengths
        self.assertEqual(1, len(measure["ALLA"][0].values))
        self.assertAlmostEqual(0.2225376, measure["ALLA"][0].values[0])
        # Second six are angles
        self.assertEqual(1, len(measure["ALLA"][6].values))
        self.assertAlmostEqual(0.2225376, measure["ALLA"][6].values[0])
        # Final six are dihedrals
        self.assertEqual(1, len(measure["ALLA"][12].values))
        self.assertAlmostEqual(0.2225376, measure["ALLA"][12].values[0])


class AtomTest(unittest.TestCase):
    def test_atom_create(self):
        atom = Atom(name="Name", num=0, typ="Type")
        self.assertEqual("Name", atom.name)
        self.assertEqual(0, atom.num)
        self.assertEqual("Type", atom.typ)


class ResidueTest(unittest.TestCase):
    def test_residue_create(self):
        residue = Residue(name="Resname")
        self.assertEqual("Resname", residue.name)

    def test_residue_add_atoms(self):
        atom = Atom(name="Name", num=0, typ="Type")
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
        frame.output_gro("water-out.gro")
        self.assertTrue(filecmp.cmp("test/data/water.gro", "water-out.gro"))

    def test_frame_read_xtc(self):
        frame = Frame(gro="test/data/water.gro", xtc="test/data/water.xtc")
        np.testing.assert_allclose(np.array([0.696, 1.33, 1.211]),
                                   frame.residues[0].atoms[0].coords)
        frame.next_frame()
        # TODO find out why this is different by 0.001
        np.testing.assert_allclose(np.array([0.696, 1.33, 1.21]),
                                   frame.residues[0].atoms[0].coords)
        frame.next_frame()
        np.testing.assert_allclose(np.array([1.176, 1.152, 1.586]),
                                   frame.residues[0].atoms[0].coords)


if __name__ == '__main__':
    unittest.main()
