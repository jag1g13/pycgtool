import unittest
import filecmp
import os
import logging

import numpy as np

from pycgtool.frame import Atom, Residue
from pycgtool.frame import NonMatchingSystemError, UnsupportedFormatException
from pycgtool.frame import Trajectory as Frame

class AtomTest(unittest.TestCase):
    def test_atom_create(self):
        atom = Atom(name="Name", num=0, type="Type")
        self.assertEqual("Name", atom.name)
        self.assertEqual(0, atom.num)
        self.assertEqual("Type", atom.type)

    def test_atom_add_missing_data(self):
        atom1 = Atom("Name1", 0, type="Type")
        atom2 = Atom("Name2", 0, mass=1)

        with self.assertRaises(AssertionError):
            atom1.add_missing_data(atom2)

        atom2 = Atom("Name1", 0, mass=1)
        atom1.add_missing_data(atom2)
        self.assertEqual(1, atom1.mass)


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


# TODO add ITP parsing tests
class FrameTest(unittest.TestCase):
    def helper_read_xtc(self, frame, first_only=False, skip_names=False):
        self.assertEqual(663, frame.natoms)

        residue = frame.residue(0)

        self.assertEqual(221, len(list(frame.residues)))
        self.assertEqual(3, residue.n_atoms)

        if not skip_names:  # MDTraj renames water
            self.assertEqual("SOL", residue.name)
            self.assertEqual("OW", residue.atom(0).name)

        atom0_coords = np.array([
            [0.696, 1.330, 1.211],
            [1.176, 1.152, 1.586],
            [1.122, 1.130, 1.534]
        ])

        box_vectors = np.array([
            [1.89868,    1.89868,    1.89868],
            [1.9052,     1.9052,     1.9052],
            [1.90325272, 1.90325272, 1.90325272]
        ])

        for i in range(1 if first_only else len(atom0_coords)):
            np.testing.assert_allclose(atom0_coords[i], residue.atom(0).coords)
            np.testing.assert_allclose(box_vectors[i], frame.box, rtol=1e-4)  # PDB files are f9.3
            frame.next_frame()

    def test_frame_create(self):
        Frame()

    def test_frame_add_residue(self):
        frame = Frame()
        residue = frame.add_residue('TEST_RESIDUE')

        self.assertTrue(residue is frame.residue(0))

    def test_frame_read_gro(self):
        # logging.disable(logging.WARNING)
        frame = Frame('test/data/water.gro')
        # logging.disable(logging.NOTSET)

        self.helper_read_xtc(frame, first_only=True, skip_names=True)

    def test_frame_read_pdb(self):
        # logging.disable(logging.WARNING)
        frame = Frame("test/data/water.pdb")
        # logging.disable(logging.NOTSET)

        self.helper_read_xtc(frame, first_only=True, skip_names=True)

    def test_frame_any_read_unsupported(self):
        with self.assertRaises(UnsupportedFormatException):
            _ = Frame('test/data/dppc.map')

    def test_frame_output_gro(self):
        frame = Frame("test/data/water.gro")
        frame.output("water-out.gro", format="gro")
        self.assertTrue(filecmp.cmp("test/data/water.gro", "water-out.gro"))
        os.remove("water-out.gro")

    def test_frame_read_xtc_numframes(self):
        logging.disable(logging.WARNING)
        frame = Frame('test/data/water.gro',
                      'test/data/water.xtc')
        logging.disable(logging.NOTSET)
        self.assertEqual(11, frame.numframes)

    def test_frame_read_xtc(self):
        logging.disable(logging.WARNING)
        frame = Frame('test/data/water.gro',
                      'test/data/water.xtc')
        logging.disable(logging.NOTSET)

        self.helper_read_xtc(frame, skip_names=True)

    def test_frame_write_xtc(self):
        try:
            os.remove("water_test2.xtc")
        except IOError:
            pass

        logging.disable(logging.WARNING)
        frame = Frame('test/data/water.gro', 'test/data/water.xtc')
        logging.disable(logging.NOTSET)

        while frame.next_frame():
            frame.write_xtc("water_test2.xtc")

    def test_raise_nonmatching_system_mdtraj(self):
        with self.assertRaises(NonMatchingSystemError):
            _ = Frame('test/data/water.gro', 'test/data/sugar.xtc')



if __name__ == '__main__':
    unittest.main()
