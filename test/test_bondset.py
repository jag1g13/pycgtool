import unittest

from pycgtool.bondset import BondSet
from pycgtool.frame import Frame
from pycgtool.mapping import Mapping


class DummyOptions:
    constr_threshold = 100000
    map_center = "geom"


class BondSetTest(unittest.TestCase):
    def test_bondset_create(self):
        measure = BondSet("test/data/sugar.bnds", DummyOptions)
        self.assertEqual(2, len(measure))
        self.assertTrue("SOL" in measure)
        self.assertTrue("ALLA" in measure)
        self.assertEqual(0, len(measure["SOL"]))
        self.assertEqual(18, len(measure["ALLA"]))

    def test_bondset_apply(self):
        measure = BondSet("test/data/sugar.bnds", DummyOptions)
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

    def test_bondset_boltzmann_invert(self):
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

        measure = BondSet("test/data/sugar.bnds", DummyOptions)
        frame = Frame("test/data/sugar.gro", xtc="test/data/sugar.xtc")
        mapping = Mapping("test/data/sugar.map", DummyOptions)
        cgframe = mapping.apply(frame)
        measure.apply(cgframe)
        frame.next_frame()
        cgframe = mapping.apply(frame)
        measure.apply(cgframe)
        measure.boltzmann_invert()
        for i, bond in enumerate(measure["ALLA"]):
            self.assertAlmostEqual(ref[i][0], bond.eqm, delta=abs(ref[i][0] / 1000000))
            self.assertAlmostEqual(ref[i][1], bond.fconst, delta=abs(ref[i][1] / 1000000))

    def test_bondset_polymer(self):
        bondset = BondSet("test/data/polyethene.bnd", DummyOptions)
        frame = Frame("test/data/polyethene.gro")
        bondset.apply(frame)
        self.assertEqual(5, len(bondset["ETH"][0].values))
        self.assertEqual(4, len(bondset["ETH"][1].values))
        self.assertEqual(4, len(bondset["ETH"][2].values))
        self.assertEqual(4, len(bondset["ETH"][3].values))
        bondset.boltzmann_invert()
        self.assertAlmostEqual(0.107, bondset["ETH"][0].eqm, places=3)
        self.assertAlmostEqual(0.107, bondset["ETH"][1].eqm, places=3)
