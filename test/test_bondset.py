import unittest

import math
import os
import pathlib
import tempfile

from pycgtool.bondset import BondSet
from pycgtool.frame import Frame
from pycgtool.mapping import Mapping
from pycgtool.util import cmp_file_whitespace_float


class DummyOptions:
    constr_threshold = 100000
    map_center = 'geom'
    virtual_map_center = 'geom'
    angle_default_fc = False
    generate_angles = True
    generate_dihedrals = False


class DummyBond:
    def __init__(self, atoms, eqm, fconst, values=None):
        self.atoms = atoms
        self.eqm = eqm
        self.fconst = fconst
        self.values = [] if values is None else values

    def __iter__(self):
        return iter(self.atoms)


# TODO add setup and teardown to put all files in a tmp directory
class BondSetTest(unittest.TestCase):
    base_dir = pathlib.Path(__file__).absolute().parent
    data_dir = base_dir.joinpath('data')

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
        ( 1.695540843207,    253.75691,   50),
        (-2.074156183261,    125.04591,   50),
        ( 2.768063749729,    135.50927,   50),
        (-2.213755550432,     51.13975,   50),
        ( 1.456495664734,     59.38162,   50),
        (-1.826134298998,    279.80889,   50)
    ]

    def test_bondset_create(self):
        measure = BondSet('test/data/sugar.bnd', DummyOptions)
        self.assertEqual(1, len(measure))
        self.assertTrue('ALLA' in measure)
        self.assertEqual(18, len(measure['ALLA']))

    def test_bondset_apply(self):
        measure = BondSet('test/data/sugar.bnd', DummyOptions)
        frame = Frame('test/data/sugar-cg.gro')
        measure.apply(frame)

        # First six are bond lengths
        self.assertEqual(1, len(measure['ALLA'][0].values))
        self.assertAlmostEqual(0.2225376, measure['ALLA'][0].values[0],
                               delta=0.2225376 / 500)
        # Second six are angles
        self.assertEqual(1, len(measure['ALLA'][6].values))
        expected = math.radians(77.22779289)
        self.assertAlmostEqual(expected, measure['ALLA'][6].values[0],
                               delta=expected / 500)
        # Final six are dihedrals
        self.assertEqual(1, len(measure['ALLA'][12].values))
        expected = math.radians(-89.5552903)
        self.assertAlmostEqual(expected, measure['ALLA'][12].values[0],
                               delta=abs(expected) / 500)

    def test_bondset_remove_triangles(self):
        bondset = BondSet('test/data/triangle.bnd', DummyOptions)
        angles = bondset.get_bond_angles('TRI', exclude_triangle=False)
        self.assertEqual(3, len(angles))
        angles = bondset.get_bond_angles('TRI', exclude_triangle=True)
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
        measure = BondSet('test/data/sugar.bnd', DummyOptions)
        frame = Frame('test/data/sugar.gro', 'test/data/sugar.xtc')
        mapping = Mapping('test/data/sugar.map', DummyOptions)

        cg_frame = mapping.apply(frame)

        measure.apply(cg_frame)
        measure.boltzmann_invert()

        self.support_check_mean_fc(measure['ALLA'], 1)

    def test_bondset_boltzmann_invert_default_fc(self):
        class DefaultOptions(DummyOptions):
            default_fc = True

        measure = BondSet('test/data/sugar.bnd', DefaultOptions)
        frame = Frame('test/data/sugar.gro', 'test/data/sugar.xtc')
        mapping = Mapping('test/data/sugar.map', DefaultOptions)

        cg_frame = mapping.apply(frame)

        measure.apply(cg_frame)
        measure.boltzmann_invert()

        self.support_check_mean_fc(measure['ALLA'], 2)

    def test_bondset_boltzmann_invert_manual_default_fc(self):
        class FuncFormOptions(DummyOptions):
            length_form = 'MartiniDefaultLength'
            angle_form = 'MartiniDefaultAngle'
            dihedral_form = 'MartiniDefaultDihedral'

        measure = BondSet('test/data/sugar.bnd', FuncFormOptions)
        frame = Frame('test/data/sugar.gro', 'test/data/sugar.xtc')
        mapping = Mapping('test/data/sugar.map', DummyOptions)

        cg_frame = mapping.apply(frame)

        measure.apply(cg_frame)
        measure.boltzmann_invert()

        self.support_check_mean_fc(measure['ALLA'], 2)

    def test_bondset_polymer(self):
        bondset = BondSet('test/data/polyethene.bnd', DummyOptions)
        frame = Frame('test/data/polyethene.gro')
        bondset.apply(frame)

        self.assertEqual(5, len(bondset['ETH'][0].values))
        self.assertEqual(4, len(bondset['ETH'][1].values))
        self.assertEqual(4, len(bondset['ETH'][2].values))
        self.assertEqual(4, len(bondset['ETH'][3].values))

        bondset.boltzmann_invert()

        self.assertAlmostEqual(0.107, bondset['ETH'][0].eqm,
                               delta=0.107 / 500)
        self.assertAlmostEqual(0.107, bondset['ETH'][1].eqm,
                               delta=0.107 / 500)

    def test_bondset_pbc(self):
        bondset = BondSet('test/data/polyethene.bnd', DummyOptions)
        frame = Frame('test/data/pbcpolyethene.gro')

        bondset.apply(frame)
        bondset.boltzmann_invert()

        for bond in bondset.get_bond_lengths('ETH', True):
            self.assertAlmostEqual(1., bond.eqm)
            self.assertEqual(float('inf'), bond.fconst)

    def test_full_itp_sugar(self):
        measure = BondSet('test/data/sugar.bnd', DummyOptions)
        frame = Frame('test/data/sugar.gro', 'test/data/sugar.xtc')
        mapping = Mapping('test/data/sugar.map', DummyOptions)

        cg_frame = mapping.apply(frame)

        measure.apply(cg_frame)
        measure.boltzmann_invert()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)

            measure.write_itp(tmp_path.joinpath('sugar_out.itp'), mapping)

            self.assertTrue(
                cmp_file_whitespace_float(
                    tmp_path.joinpath('sugar_out.itp'),
                    self.data_dir.joinpath('sugar_out.itp'),
                    rtol=0.005,
                    verbose=True))

    def test_full_itp_vsites(self):
        """Test full operation to output of .itp file for molecule with vsites."""
        options = DummyOptions()
        options.generate_angles = False

        measure = BondSet('test/data/martini3/naphthalene.bnd', options)
        frame = Frame('test/data/martini3/naphthalene.gro')
        mapping = Mapping('test/data/martini3/naphthalene.map', options)

        cg_frame = mapping.apply(frame)

        measure.apply(cg_frame)
        measure.boltzmann_invert()

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            measure.write_itp(tmp_path.joinpath('naphthalene_out.itp'),
                              mapping)

            self.assertTrue(
                cmp_file_whitespace_float(
                    tmp_path.joinpath('naphthalene_out.itp'),
                    self.data_dir.joinpath('martini3/naphthalene_out.itp'),
                    rtol=0.005,
                    verbose=True))

    def test_duplicate_atoms_in_bond(self):
        with self.assertRaises(ValueError):
            _ = BondSet('test/data/duplicate_atoms.bnd', DummyOptions)

    def test_dump_bonds(self):
        measure = BondSet('test/data/sugar.bnd', DummyOptions)
        frame = Frame('test/data/sugar.gro', 'test/data/sugar.xtc')
        mapping = Mapping('test/data/sugar.map', DummyOptions)

        cg_frame = mapping.apply(frame)

        measure.apply(cg_frame)
        measure.boltzmann_invert()
        measure.dump_values()

        data_dir = pathlib.Path('test/data')
        filenames = ('ALLA_length.dat', 'ALLA_angle.dat', 'ALLA_dihedral.dat')
        for filename in filenames:
            self.assertTrue(
                cmp_file_whitespace_float(data_dir.joinpath(filename),
                                          filename,
                                          rtol=0.008,
                                          verbose=True))
            os.remove(filename)

    def test_get_lines_for_bond_dump(self):
        expected = [
            '     0.00000     1.00000     2.00000',
            '     1.00000     2.00000     3.00000',
            '     2.00000     3.00000     4.00000',
            '     3.00000     4.00000     5.00000'
        ]

        bonds = [
            DummyBond(None, None, None, values=[0, 1, 2, 3]),
            DummyBond(None, None, None, values=[1, 2, 3, 4]),
            DummyBond(None, None, None, values=[2, 3, 4, 5])
        ]

        output = BondSet._get_lines_for_bond_dump(bonds)

        self.assertListEqual(expected, output)

    def test_get_lines_for_bond_dump_sample(self):
        expected = [
            '     0.00000     1.00000     2.00000',
            '     1.00000     2.00000     3.00000',
            '     2.00000     3.00000     4.00000',
            '     3.00000     4.00000     5.00000'
        ]

        bonds = [
            DummyBond(None, None, None, values=[0, 1, 2, 3]),
            DummyBond(None, None, None, values=[1, 2, 3, 4]),
            DummyBond(None, None, None, values=[2, 3, 4, 5])
        ]

        nlines = 2
        output = BondSet._get_lines_for_bond_dump(bonds, target_number=nlines)

        self.assertEqual(nlines, len(output))

        seen = set()
        for line in output:
            self.assertIn(line, expected)
            self.assertNotIn(line, seen)
            seen.add(line)
