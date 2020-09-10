import unittest
import filecmp
import os

import numpy as np

from pycgtool.mapping import Mapping, VirtualMap
from pycgtool.frame import Frame


class DummyOptions:
    map_center = "geom"
    virtual_map_center = "geom"


class MappingTest(unittest.TestCase):
    def test_mapping_create(self):
        """Test that a mapping can be correctly read from a file."""
        mapping = Mapping('test/data/water.map', DummyOptions)
        self.assertEqual(2, len(mapping))  # SOL and HOH
        self.assertTrue('SOL' in mapping)

        sol_map = mapping['SOL']

        self.assertEqual(1, len(sol_map))
        self.assertEqual(3, len(sol_map[0].atoms))
        self.assertEqual('OW', sol_map[0].atoms[0])
        self.assertEqual('HW1', sol_map[0].atoms[1])
        self.assertEqual('HW2', sol_map[0].atoms[2])

    def test_mapping_rename(self):
        """Test that an alternate mapping is created using MDTraj conventions."""
        mapping = Mapping('test/data/water.map', DummyOptions)
        self.assertEqual(2, len(mapping))  # SOL and HOH
        self.assertTrue('HOH' in mapping)

        sol_map = mapping['HOH']

        self.assertEqual(1, len(sol_map))
        self.assertEqual(3, len(sol_map[0].atoms))
        self.assertEqual('O', sol_map[0].atoms[0])
        self.assertEqual('H1', sol_map[0].atoms[1])
        self.assertEqual('H2', sol_map[0].atoms[2])

    def test_virtual_mapping_create(self):
        mapping = Mapping("test/data/martini3/naphthalene.map", DummyOptions)
        self.assertEqual(1, len(mapping))
        self.assertTrue("NAPH" in mapping)
        self.assertEqual(5, len(mapping["NAPH"]))
        self.assertTrue(isinstance(mapping["NAPH"][2], VirtualMap))
        self.assertEqual(4, len(mapping["NAPH"][2].atoms))
        self.assertEqual("R1", mapping["NAPH"][2].atoms[0])
        self.assertEqual(1, [isinstance(bead, VirtualMap) for bead in mapping["NAPH"]].count(True))

    def test_mapping_apply(self):
        mapping = Mapping("test/data/water.map", DummyOptions)
        frame = Frame("test/data/water.gro")
        cg_frame = mapping.apply(frame)

        self.assertEqual(frame.natoms / 3, cg_frame.natoms)

        cg_frame.save("water-cg.gro")

        self.assertTrue(filecmp.cmp("test/data/water-cg.gro", "water-cg.gro"))
        os.remove("water-cg.gro")

    def test_mapping_charges(self):
        mapping = Mapping("test/data/dppc.map", DummyOptions)
        self.assertEqual( 1, mapping["DPPC"][0].charge)
        self.assertEqual(-1, mapping["DPPC"][1].charge)

    def test_mapping_pbc(self):
        frame = Frame("test/data/pbcwater.gro")

        mapping = Mapping("test/data/water.map", DummyOptions)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(frame.atom(0).coords, cg_frame.atom(0).coords)

    def test_mapping_weights_geom(self):
        frame = Frame("test/data/two.gro")

        mapping = Mapping("test/data/two.map", DummyOptions)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(np.array([1.5, 1.5, 1.5]),
                                   cg_frame.residue(0).atom(0).coords)

    def test_virtual_mapping_weights_geom(self):
        frame = Frame("test/data/martini3/four.gro")

        mapping = Mapping("test/data/martini3/four.map", DummyOptions)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(np.array([2.5, 2.5, 2.5]), 
                                   cg_frame.residue(0).atom(2).coords)

    def test_mapping_weights_mass(self):
        frame = Frame("test/data/two.gro")
        options = DummyOptions()
        options.map_center = "mass"

        mapping = Mapping("test/data/two.map", options, itp_filename="test/data/two.itp")
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(np.array([2., 2., 2.]),
                                   cg_frame.residue(0).atom(0).coords)

    def test_virtual_mapping_weights_mass(self):
        frame = Frame("test/data/martini3/four.gro")
        options = DummyOptions()
        options.virtual_map_center = "mass"

        mapping = Mapping("test/data/martini3/four.map", options, itp_filename="test/data/martini3/four.itp")
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(np.array([3.0, 3.0, 3.0]),
                                   cg_frame.residue(0).atom(2).coords)

    def test_mapping_weights_guess_mass(self):
        frame = Frame("test/data/two.gro")
        options = DummyOptions()
        options.map_center = "mass"

        mapping = Mapping("test/data/two.map", options)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(np.array([1.922575,  1.922575,  1.922575], dtype=np.float32),
                                   cg_frame.residue(0).atom(0).coords, rtol=0.01)

    def test_virtual_mapping_weights_guess_mass(self):
        frame = Frame("test/data/martini3/four.gro")
        options = DummyOptions()
        options.virtual_map_center = "mass"

        mapping = Mapping("test/data/martini3/four.map", options)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(np.array([2.83337, 2.83337, 2.83337], dtype=np.float32),
                                   cg_frame.residue(0).atom(2).coords, rtol=0.01)

    def test_mapping_weights_first(self):
        frame = Frame("test/data/two.gro")
        options = DummyOptions()
        options.map_center = "first"

        mapping = Mapping("test/data/two.map", options, itp_filename="test/data/two.itp")
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(np.array([1., 1., 1.]),
                                   cg_frame.residue(0).atom(0).coords)

    def test_mapping_itp_multi(self):
        mapping = Mapping("test/data/membrane/membrane.map",
                          DummyOptions,
                          itp_filename="test/data/membrane/membrane.top")
        self.assertAlmostEqual( -1.2, mapping["POPE"][0].charge, delta=0.0001)
        self.assertAlmostEqual(0, mapping["POPG"][0].charge, delta=0.0001)

        self.assertAlmostEqual([94.9716], mapping["POPE"][0].mass, delta=0.0001)
