import unittest

from pycgtool.bondset import BondSet
from pycgtool.frame import Frame
from pycgtool.mapping import Mapping
from pycgtool.util import cmp_whitespace_float


class DummyOptions:
    constr_threshold = 100000
    map_center = "geom"
    angle_default_fc = False
    generate_angles = True
    generate_dihedrals = False


class BondSetTest(unittest.TestCase):
    def test_bondset_create(self):
        measure = BondSet("test/data/sugar.bnd", DummyOptions)
        self.assertEqual(2, len(measure))
        self.assertTrue("SOL" in measure)
        self.assertTrue("ALLA" in measure)
        self.assertEqual(0, len(measure["SOL"]))
        self.assertEqual(18, len(measure["ALLA"]))

    def test_bondset_apply(self):
        measure = BondSet("test/data/sugar.bnd", DummyOptions)
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
        ref = [(  0.222817161647,   9560719),
               (  0.216766804211,  27930081),
               (  0.223815134299,   3099614),
               (  0.239439104442, 322294637),
               (  0.215497208685,   3683331),
               (  0.179544192953,  10260492),
               ( 76.9617550012,    13048242),
               (116.009786726,       927619),
               (109.247878273,       195899),
               ( 84.4955293511,      871991),
               (149.557001531,       471274),
               ( 99.0524708671,      100943),
               (-89.5321545897,       91248),
               ( 70.040930903,       169002),
               (-21.1194544981,      995309),
               ( 48.3727747848,       35712),
               (-85.923133935,        62394),
               ( 70.3195444564,      778389)]

        measure = BondSet("test/data/sugar.bnd", DummyOptions)
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
            # Allows for slight modifications to code
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
            self.assertEqual(float("inf"), bond.fconst)

    def test_full_itp_sugar(self):
        measure = BondSet("test/data/sugar.bnd", DummyOptions)
        frame = Frame("test/data/sugar.gro", xtc="test/data/sugar.xtc")
        mapping = Mapping("test/data/sugar.map", DummyOptions)
        cgframe = mapping.apply(frame)

        while frame.next_frame():
            cgframe = mapping.apply(frame, cgframe=cgframe, exclude={"SOL"})
            measure.apply(cgframe)

        measure.boltzmann_invert()
        measure.write_itp("sugar_out.itp", mapping, exclude={"SOL"})

        self.assertTrue(cmp_whitespace_float("sugar_out.itp", "test/data/sugar_out.itp", float_rel_error=0.001))
