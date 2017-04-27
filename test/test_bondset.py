import unittest

import logging
import math

from pycgtool.bondset import BondSet
from pycgtool.frame import Frame
from pycgtool.mapping import Mapping
from pycgtool.util import cmp_whitespace_float

try:
    import mdtraj
    mdtraj_present = True
except ImportError:
    mdtraj_present = False


class DummyOptions:
    constr_threshold = 100000
    map_center = "geom"
    angle_default_fc = False
    generate_angles = True
    generate_dihedrals = False


class BondSetTest(unittest.TestCase):
    # Columns are: eqm value, standard fc, MARTINI defaults fc
    invert_test_ref_data = [
        ( 0.220474419132,  72658.18163, 1250),
        ( 0.221844516963,  64300.01188, 1250),
        ( 0.216313356504,  67934.93368, 1250),
        ( 0.253166204438,  19545.27388, 1250),
        ( 0.205958461836,  55359.06367, 1250),
        ( 0.180550961226, 139643.66601, 1250),
        ( 1.359314249473,    503.24211,   25),
        ( 2.026002746003,    837.76904,   25),
        ( 1.937848056214,    732.87969,   25),
        ( 1.453592079716,    945.32633,   25),
        ( 2.504189933347,    771.63691,   25),
        ( 1.733002945619,    799.82825,   25),
        (-1.446051810383,    253.75691,   50),
        ( 1.067436470329,    125.04591,   50),
        (-0.373528903861,    135.50927,   50),
        ( 0.927837103158,     51.13975,   50),
        (-1.685096988856,     59.38162,   50),
        ( 1.315458354592,    279.80889,   50)
    ]

    def test_bondset_create(self):
        measure = BondSet("test/data/sugar.bnd", DummyOptions)
        self.assertEqual(1, len(measure))
        self.assertTrue("ALLA" in measure)
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
        expected = math.radians(77.22779289)
        self.assertAlmostEqual(expected, measure["ALLA"][6].values[0],
                               delta=expected / 500)
        # Final six are dihedrals
        self.assertEqual(1, len(measure["ALLA"][12].values))
        expected = math.radians(-89.5552903)
        self.assertAlmostEqual(expected, measure["ALLA"][12].values[0],
                               delta=abs(expected) / 500)

    def test_bondset_remove_triangles(self):
        bondset = BondSet("test/data/triangle.bnd", DummyOptions)
        angles = bondset.get_bond_angles("TRI", exclude_triangle=False)
        self.assertEqual(3, len(angles))
        angles = bondset.get_bond_angles("TRI", exclude_triangle=True)
        self.assertEqual(0, len(angles))

    def support_check_mean_fc(self, mol_bonds, fc_column_number):
        # Require accuracy to 0.5%
        # Allows for slight modifications to code
        accuracy = 0.005

        for i, bond in enumerate(mol_bonds):
            ref = self.invert_test_ref_data
            self.assertAlmostEqual(ref[i][0], bond.eqm,
                                   delta=abs(ref[i][0] * accuracy))
            self.assertAlmostEqual(ref[i][fc_column_number], bond.fconst,
                                   delta=abs(ref[i][fc_column_number] * accuracy))

    def test_bondset_boltzmann_invert(self):
        measure = BondSet("test/data/sugar.bnd", DummyOptions)
        frame = Frame("test/data/sugar.gro", xtc="test/data/sugar.xtc")
        mapping = Mapping("test/data/sugar.map", DummyOptions)

        cgframe = mapping.apply(frame)
        while frame.next_frame():
            cgframe = mapping.apply(frame, cgframe=cgframe)
            measure.apply(cgframe)

        measure.boltzmann_invert()
        self.support_check_mean_fc(measure["ALLA"], 1)

    def test_bondset_boltzmann_invert_default_fc(self):
        class DefaultOptions(DummyOptions):
            default_fc = True

        measure = BondSet("test/data/sugar.bnd", DefaultOptions)
        frame = Frame("test/data/sugar.gro", xtc="test/data/sugar.xtc")
        mapping = Mapping("test/data/sugar.map", DefaultOptions)

        cgframe = mapping.apply(frame)
        while frame.next_frame():
            cgframe = mapping.apply(frame, cgframe=cgframe)
            measure.apply(cgframe)

        measure.boltzmann_invert()
        self.support_check_mean_fc(measure["ALLA"], 2)

    def test_bondset_boltzmann_invert_manual_default_fc(self):
        class FuncFormOptions(DummyOptions):
            length_form = "MartiniDefaultLength"
            angle_form = "MartiniDefaultAngle"
            dihedral_form = "MartiniDefaultDihedral"

        measure = BondSet("test/data/sugar.bnd", FuncFormOptions)
        frame = Frame("test/data/sugar.gro", xtc="test/data/sugar.xtc")
        mapping = Mapping("test/data/sugar.map", DummyOptions)

        cgframe = mapping.apply(frame)
        while frame.next_frame():
            cgframe = mapping.apply(frame, cgframe=cgframe)
            measure.apply(cgframe)

        measure.boltzmann_invert()
        self.support_check_mean_fc(measure["ALLA"], 2)

    @unittest.skipIf(not mdtraj_present, "MDTRAJ or Scipy not present")
    def test_bondset_boltzmann_invert_mdtraj(self):
        logging.disable(logging.WARNING)
        frame = Frame("test/data/sugar.gro", xtc="test/data/sugar.xtc",
                      xtc_reader="mdtraj")
        logging.disable(logging.NOTSET)

        measure = BondSet("test/data/sugar.bnd", DummyOptions)
        mapping = Mapping("test/data/sugar.map", DummyOptions)

        cgframe = mapping.apply(frame)
        while frame.next_frame():
            cgframe = mapping.apply(frame, cgframe=cgframe)
            measure.apply(cgframe)

        measure.boltzmann_invert()
        self.support_check_mean_fc(measure["ALLA"], 1)

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
            cgframe = mapping.apply(frame, cgframe=cgframe)
            measure.apply(cgframe)

        measure.boltzmann_invert()

        logging.disable(logging.WARNING)
        measure.write_itp("sugar_out.itp", mapping)
        logging.disable(logging.NOTSET)

        self.assertTrue(cmp_whitespace_float("sugar_out.itp", "test/data/sugar_out.itp", float_rel_error=0.001))

    def test_duplicate_atoms_in_bond(self):
        with self.assertRaises(ValueError):
            bondset = BondSet("test/data/duplicate_atoms.bnd", DummyOptions)
