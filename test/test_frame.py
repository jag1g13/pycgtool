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

    def test_measure_apply(self):
        measure = Measure("test/data/sugar.bnds")
        frame = Frame("test/data/sugar-cg.gro")
        measure.apply(frame)
        # First six are bond lengths
        self.assertEqual(1, len(measure["ALLA"][0].values))
        self.assertAlmostEqual(0.2225376, measure["ALLA"][0].values[0])
        # Second six are angles
        self.assertEqual(1, len(measure["ALLA"][6].values))
        self.assertAlmostEqual(77.22779289, measure["ALLA"][6].values[0])
        # Final six are dihedrals
        self.assertEqual(1, len(measure["ALLA"][12].values))
        self.assertAlmostEqual(-89.5552903, measure["ALLA"][12].values[0])

    def test_measure_boltzmann_invert(self):
        ref = [(0.222817161647, 19116816.3363),
               (0.216766804211, 55835643.1619),
               (0.223815134299, 6199518.81029),
               (0.239439149124, 643501628.698),
               (0.215497173862, 7365459.20621),
               (0.179544192953, 20521533.3397),
               (76.9617550012, 26029820.9708),
               (116.009786726, 1854775.98493),
               (109.247869242, 391862.51047),
               (84.4955258944, 1744922.00401),
               (149.557009586, 942386.47818),
               (99.0524708671, 201901.375283),
               (-89.5321545897, 182481.962675),
               (70.04092798, 337927.889101),
               (-21.1194542129, 1990800.34058),
               (48.3727926542, 71431.9545811),
               (-85.9231539259, 124785.032295),
               (70.3195444564, 1555761.24713)]

        measure = Measure("test/data/sugar.bnds")
        frame = Frame("test/data/sugar.gro", xtc="test/data/sugar.xtc")
        mapping = Mapping("test/data/sugar.map")
        cgframe = mapping.apply(frame)
        measure.apply(cgframe)
        frame.next_frame()
        cgframe = mapping.apply(frame)
        measure.apply(cgframe)
        measure.boltzmann_invert()
        for i, bond in enumerate(measure["ALLA"]):
            self.assertAlmostEqual(ref[i][0], bond.eqm, delta=abs(ref[i][0] / 1000000))
            self.assertAlmostEqual(ref[i][1], bond.fconst, delta=abs(ref[i][1] / 1000000))


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
        frame.output_gro("water-out.gro")
        self.assertTrue(filecmp.cmp("test/data/water.gro", "water-out.gro"))

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
