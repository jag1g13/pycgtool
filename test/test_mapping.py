import unittest
import filecmp

from pycgtool.mapping import Mapping
from pycgtool.frame import Frame


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
        cgframe.output("water-cg.gro", format="gro")
        self.assertTrue(filecmp.cmp("test/data/water-cg.gro", "water-cg.gro"))
