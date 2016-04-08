import unittest
import filecmp
import os

import numpy as np

from pycgtool.mapping import Mapping
from pycgtool.frame import Frame, Atom


class DummyOptions:
    map_center = "geom"


class DummyOptionsMass:
    map_center = "mass"


class MappingTest(unittest.TestCase):
    def test_mapping_create(self):
        mapping = Mapping("test/data/water.map", DummyOptions)
        self.assertEqual(1, len(mapping))
        self.assertTrue("SOL" in mapping)
        self.assertEqual(1, len(mapping["SOL"]))
        self.assertEqual(3, len(mapping["SOL"][0].atoms))
        self.assertEqual("OW", mapping["SOL"][0].atoms[0])

    def test_mapping_apply(self):
        mapping = Mapping("test/data/water.map", DummyOptions)
        frame = Frame("test/data/water.gro")
        cgframe = mapping.apply(frame)
        self.assertEqual(len(frame), len(cgframe))
        cgframe.output("water-cg.gro", format="gro")
        self.assertTrue(filecmp.cmp("test/data/water-cg.gro", "water-cg.gro"))
        os.remove("water-cg.gro")

    def test_mapping_pbc(self):
        mapping = Mapping("test/data/water.map", DummyOptions)
        frame = Frame("test/data/pbcwater.gro")
        cgframe = mapping.apply(frame)
        np.testing.assert_allclose(frame[0][0].coords, cgframe[0][0].coords)

    def test_mapping_weights(self):
        frame = Frame("test/data/two.gro")

        mapping = Mapping("test/data/two.map", DummyOptions, itp="test/data/two.itp")
        cg = mapping.apply(frame)
        np.testing.assert_allclose(np.array([1.5, 1.5, 1.5]), cg[0][0].coords)

        mapping = Mapping("test/data/two.map", DummyOptionsMass, itp="test/data/two.itp")
        cg = mapping.apply(frame)
        np.testing.assert_allclose(np.array([2., 2., 2.]), cg[0][0].coords)



