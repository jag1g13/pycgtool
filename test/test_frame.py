import filecmp
import os
import unittest

import numpy as np

from pycgtool.frame import NonMatchingSystemError, UnsupportedFormatException
from pycgtool.frame import Trajectory as Frame


def try_remove(filename) -> None:
    try:
        os.remove(filename)

    except IOError:
        pass


# TODO add ITP parsing tests
class FrameTest(unittest.TestCase):
    def check_reference_topology(self, frame, skip_names=True):
        self.assertEqual(663, frame.natoms)

        residue = frame.residue(0)

        self.assertEqual(221, len(list(frame.residues)))
        self.assertEqual(3, residue.n_atoms)

        if not skip_names:  # MDTraj renames water
            self.assertEqual('SOL', residue.name)
            self.assertEqual('OW', residue.atom(0).name)

    def check_reference_frame(self, frame):
        self.check_reference_topology(frame)

        atom0_coords = np.array([0.696, 1.330, 1.211])
        box_vectors = np.array([1.89868, 1.89868, 1.89868])

        np.testing.assert_allclose(atom0_coords, frame.atom(0).coords)
        np.testing.assert_allclose(box_vectors, frame.box,
                                   rtol=1e-4)  # PDB files are f9.3

    def check_reference_trajectory(self, frame):
        self.check_reference_topology(frame)

        atom0_coords_array = np.array([[1.176, 1.152, 1.586],
                                       [1.122, 1.130, 1.534]])

        box_vectors_array = np.array([[1.9052, 1.9052, 1.9052],
                                      [1.90325272, 1.90325272, 1.90325272]])

        for _, (atom0_coords, box_vectors) in enumerate(
                zip(atom0_coords_array, box_vectors_array)):
            np.testing.assert_allclose(atom0_coords, frame.atom(0).coords)
            np.testing.assert_allclose(box_vectors, frame.box,
                                       rtol=1e-4)  # PDB files are f9.3
            frame.next_frame()

    def test_frame_create(self):
        Frame()

    def test_frame_add_residue(self):
        frame = Frame()
        residue = frame.add_residue('TEST_RESIDUE')

        self.assertTrue(residue is frame.residue(0))

    def test_frame_read_gro(self):
        frame = Frame('test/data/water.gro')

        self.check_reference_frame(frame)

    def test_frame_read_pdb(self):
        frame = Frame('test/data/water.pdb')

        self.check_reference_frame(frame)

    def test_frame_any_read_unsupported(self):
        with self.assertRaises(UnsupportedFormatException):
            _ = Frame('test/data/dppc.map')

    def test_frame_output_gro(self):
        frame = Frame('test/data/water.gro')
        frame.save('water-out.gro')
        self.assertTrue(filecmp.cmp('test/data/water.gro', 'water-out.gro'))

        try_remove('water-out.gro')

    def test_frame_read_xtc_numframes(self):
        frame = Frame('test/data/water.gro', 'test/data/water.xtc')
        self.assertEqual(10, frame.numframes)

    def test_frame_read_xtc(self):
        frame = Frame('test/data/water.gro', 'test/data/water.xtc')

        self.check_reference_trajectory(frame)

    def test_frame_write_xtc(self):
        try_remove('water_test2.xtc')

        frame = Frame('test/data/water.gro', 'test/data/water.xtc')
        frame.save('water_test2.xtc')

    def test_raise_nonmatching_system_mdtraj(self):
        with self.assertRaises(NonMatchingSystemError):
            _ = Frame('test/data/water.gro', 'test/data/sugar.xtc')


if __name__ == '__main__':
    unittest.main()
