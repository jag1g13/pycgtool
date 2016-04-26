import unittest

from pycgtool.bondset import BondSet
from pycgtool.frame import Frame
from pycgtool.mapping import Mapping


class DummyOptions:
    constr_threshold = 100000
    map_center = "geom"
    angle_default_fc = False


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
        self.assertAlmostEqual(0.2225376, measure["ALLA"][0].values[0],
                               delta=0.2225376 / 500)
        # Second six are angles
        self.assertEqual(1, len(measure["ALLA"][6].values))
        self.assertAlmostEqual(77.22779289, measure["ALLA"][6].values[0],
                               delta=77.22779289 / 500)
        # Final six are dihedrals
        self.assertEqual(1, len(measure["ALLA"][12].values))
        self.assertAlmostEqual(-89.5552903, measure["ALLA"][12].values[0],
                               delta=89.552903 / 500)

    def test_bondset_remove_triangles(self):
        bondset = BondSet("test/data/triangle.bnd", DummyOptions)
        angles = bondset.get_bond_angles("TRI", exclude_triangle=False)
        self.assertEqual(3, len(angles))
        angles = bondset.get_bond_angles("TRI", exclude_triangle=True)
        self.assertEqual(0, len(angles))

    def test_bondset_boltzmann_invert(self):
        ref = [(0.222817161647, 19116816.3363),
               (0.216766804211, 55835643.1619),
               (0.223815134299, 6199518.81029),
               (0.239439104442, 642593936.491),
               (0.215497208685, 7366326.46847),
               (0.179544192953, 20521533.3397),
               (76.9617550012, 23447621.6268),
               (116.009786726, 1210005.81229),
               (109.247878273, 311284.586971),
               (84.4955293511, 1713127.90333),
               (149.557001531, 62116.1326728),
               (99.0524708671, 192028.755165),
               (-89.5321545897, 182481.962675),
               (70.040930903, 337940.374373),
               (-21.1194544981, 1990817.75509),
               (48.3727747848, 71424.537429),
               (-85.923133935, 124765.874439),
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
            # Require accuracy to 0.5%
            # Allows for slight modifications to algorithm
            self.assertAlmostEqual(ref[i][0], bond.eqm, delta=abs(ref[i][0] / 200))
            self.assertAlmostEqual(ref[i][1], bond.fconst, delta=abs(ref[i][1] / 200))

    def test_bondset_polymer(self):
        bondset = BondSet("test/data/polyethene.bnd", DummyOptions)
        frame = Frame("test/data/polyethene.gro")
        bondset.apply(frame)
        self.assertEqual(5, len(bondset["ETH"][0].values))
        self.assertEqual(4, len(bondset["ETH"][1].values))
        self.assertEqual(4, len(bondset["ETH"][2].values))
        self.assertEqual(4, len(bondset["ETH"][3].values))
        bondset.boltzmann_invert()
        self.assertAlmostEqual(0.107, bondset["ETH"][0].eqm,
                               delta=0.107 / 500)
        self.assertAlmostEqual(0.107, bondset["ETH"][1].eqm,
                               delta=0.107 / 500)

    def test_bondset_pbc(self):
        bondset = BondSet("test/data/polyethene.bnd", DummyOptions)
        frame = Frame("test/data/pbcpolyethene.gro")
        bondset.apply(frame)
        bondset.boltzmann_invert()
        for bond in bondset.get_bond_lengths("ETH", True):
            self.assertAlmostEqual(1., bond.eqm)
            self.assertEqual(0., bond.fconst)
