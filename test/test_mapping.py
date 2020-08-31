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
        mapping = Mapping("test/data/water.map", DummyOptions)
        self.assertEqual(1, len(mapping))
        self.assertTrue("SOL" in mapping)
        self.assertEqual(1, len(mapping["SOL"]))
        self.assertEqual(3, len(mapping["SOL"][0].atoms))
        self.assertEqual("OW", mapping["SOL"][0].atoms[0])

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
        cgframe = mapping.apply(frame)
        self.assertEqual(len(frame), len(cgframe))
        cgframe.output("water-cg.gro", format="gro")
        self.assertTrue(filecmp.cmp("test/data/water-cg.gro", "water-cg.gro"))
        os.remove("water-cg.gro")

    def test_mapping_charges(self):
        mapping = Mapping("test/data/dppc.map", DummyOptions)
        self.assertEqual( 1, mapping["DPPC"][0].charge)
        self.assertEqual(-1, mapping["DPPC"][1].charge)
        frame = Frame("test/data/dppc.gro")
        cgframe = mapping.apply(frame)
        self.assertEqual( 1, cgframe[0][0].charge)
        self.assertEqual(-1, cgframe[0][1].charge)

    def test_mapping_pbc(self):
        mapping = Mapping("test/data/water.map", DummyOptions)
        frame = Frame("test/data/pbcwater.gro")
        cgframe = mapping.apply(frame)
        np.testing.assert_allclose(frame[0][0].coords, cgframe[0][0].coords)

    def test_mapping_weights_geom(self):
        frame = Frame("test/data/two.gro")
        mapping = Mapping("test/data/two.map", DummyOptions)
        cg = mapping.apply(frame)
        np.testing.assert_allclose(np.array([1.5, 1.5, 1.5]), cg[0][0].coords)

    def test_virtual_mapping_weights_geom(self):
        frame = Frame("test/data/martini3/four.gro")
        mapping = Mapping("test/data/martini3/four.map", DummyOptions)
        cg = mapping.apply(frame)
        np.testing.assert_allclose(np.array([2.5, 2.5, 2.5]), cg[0][2].coords)

    def test_mapping_weights_mass(self):
        frame = Frame("test/data/two.gro")
        options = DummyOptions()
        options.map_center = "mass"

        mapping = Mapping("test/data/two.map", options, itp="test/data/two.itp")
        cg = mapping.apply(frame)
        np.testing.assert_allclose(np.array([2., 2., 2.]), cg[0][0].coords)

    def test_virtual_mapping_weights_mass(self):
        frame = Frame("test/data/martini3/four.gro")
        options = DummyOptions()
        options.virtual_map_center = "mass"
        mapping = Mapping("test/data/martini3/four.map", options, itp="test/data/martini3/four.itp")
        cg = mapping.apply(frame)
        np.testing.assert_allclose(np.array([3.0, 3.0, 3.0]), cg[0][2].coords)

    def test_mapping_weights_guess_mass(self):
        frame = Frame("test/data/two.gro")
        options = DummyOptions()
        options.map_center = "mass"

        mapping = Mapping("test/data/two.map", options)
        cg = mapping.apply(frame)
        np.testing.assert_allclose(np.array([1.922575,  1.922575,  1.922575], dtype=np.float32),
                                   cg[0][0].coords, rtol=0.01)

    def test_virtual_mapping_weights_guess_mass(self):
        frame = Frame("test/data/martini3/four.gro")
        options = DummyOptions()
        options.virtual_map_center = "mass"
        mapping = Mapping("test/data/martini3/four.map", options)
        cg = mapping.apply(frame)
        np.testing.assert_allclose(np.array([2.83337, 2.83337, 2.83337], dtype=np.float32),
                                   cg[0][2].coords, rtol=0.01)

    def test_mapping_weights_first(self):
        frame = Frame("test/data/two.gro")
        options = DummyOptions()
        options.map_center = "first"

        mapping = Mapping("test/data/two.map", options, itp="test/data/two.itp")
        cg = mapping.apply(frame)
        np.testing.assert_allclose(np.array([1., 1., 1.]), cg[0][0].coords)
