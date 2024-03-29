import filecmp
import os
import pathlib
import unittest
import tempfile

import numpy as np

from pycgtool.frame import NonMatchingSystemError, UnsupportedFormatException
from pycgtool.frame import Frame
from pycgtool import util


def try_remove(filename) -> None:
    try:
        os.remove(filename)

    except IOError:
        pass


# TODO add ITP parsing tests
class FrameTest(unittest.TestCase):
    base_dir = pathlib.Path(__file__).absolute().parent
    data_dir = base_dir.joinpath("data")

    def check_reference_topology(self, frame, skip_names=True):
        self.assertEqual(663, frame.natoms)

        residue = frame.residue(0)

        self.assertEqual(221, len(list(frame.residues)))
        self.assertEqual(3, residue.n_atoms)

        if not skip_names:  # MDTraj renames water
            self.assertEqual("SOL", residue.name)
            self.assertEqual("OW", residue.atom(0).name)

    def check_reference_frame(self, frame, check_box: bool = True):
        self.check_reference_topology(frame)

        atom0_coords = np.array([[0.696, 1.330, 1.211]])
        box_vectors = np.array([[1.89868, 1.89868, 1.89868]])

        np.testing.assert_allclose(atom0_coords, frame.atom(0).coords)
        if check_box:
            np.testing.assert_allclose(
                box_vectors, frame.unitcell_lengths, rtol=1e-4
            )  # PDB files are f9.3

    def check_reference_trajectory(self, frame):
        self.check_reference_topology(frame)
        atom0_coords = np.array(
            [
                [1.176     , 1.1520001 , 1.5860001 ],
                [1.1220001 , 1.13      , 1.534     ],
                [1.0580001 , 1.1620001 , 1.462     ],
                [0.91600007, 1.276     , 1.5580001 ],
                [0.73200005, 1.1240001 , 1.286     ],
                [0.64000005, 1.2160001 , 1.258     ],
                [0.632     , 1.312     , 1.2520001 ],
                [0.606     , 1.284     , 1.246     ],
                [0.582     , 1.312     , 1.1600001 ],
                [0.68200004, 1.22      , 1.25      ],
                [0.69600004, 1.33      , 1.21      ],
            ],
            dtype=np.float32,
        )  # fmt: skip
        unitcell_lengths = np.array(
            [
                [1.9052   , 1.9052   , 1.9052   ],
                [1.9032527, 1.9032527, 1.9032527],
                [1.9040661, 1.9040661, 1.9040661],
                [1.896811 , 1.896811 , 1.896811 ],
                [1.8985983, 1.8985983, 1.8985983],
                [1.9033976, 1.9033976, 1.9033976],
                [1.8904614, 1.8904614, 1.8904614],
                [1.9013108, 1.9013108, 1.9013108],
                [1.8946321, 1.8946321, 1.8946321],
                [1.898091 , 1.898091 , 1.898091 ],
                [1.898684 , 1.898684 , 1.898684 ],
            ],
            dtype=np.float32,
        )  # fmt: skip

        np.testing.assert_allclose(atom0_coords, frame.atom(0).coords)
        np.testing.assert_allclose(
            unitcell_lengths, frame.unitcell_lengths, rtol=1e-4
        )  # PDB files are f9.3

    def test_frame_create(self):
        Frame()

    def test_frame_add_residue(self):
        frame = Frame()
        residue = frame.add_residue("TEST_RESIDUE")

        self.assertTrue(residue is frame.residue(0))

    def test_frame_read_gro(self):
        frame = Frame(self.data_dir.joinpath("water.gro"))

        self.check_reference_frame(frame)

    def test_frame_read_pdb(self):
        frame = Frame(self.data_dir.joinpath("water.pdb"))

        self.check_reference_frame(frame)

    def test_frame_read_pdb_nobox(self):
        """Test reading a PDB with zero box volume."""
        frame = Frame(self.data_dir.joinpath("water-nobox.pdb"))

        self.check_reference_frame(frame, check_box=False)

    def test_frame_read_zero_box(self):
        frame = Frame(self.data_dir.joinpath("polyethene.gro"))
        self.assertIsNone(frame.unitcell_lengths)

    def test_frame_any_read_unsupported(self):
        with self.assertRaises(UnsupportedFormatException):
            _ = Frame(self.data_dir.joinpath("dppc.map"))

    def test_frame_output_gro(self):
        frame = Frame(self.data_dir.joinpath("water.gro"))

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)

            frame.save(tmp_path.joinpath("water-out.gro"))
            self.assertTrue(
                filecmp.cmp(
                    self.data_dir.joinpath("water.gro"),
                    tmp_path.joinpath("water-out.gro"),
                )
            )

    def test_frame_read_xtc_numframes(self):
        frame = Frame(
            self.data_dir.joinpath("water.gro"), self.data_dir.joinpath("water.xtc")
        )
        self.assertEqual(11, frame.n_frames)

    def test_frame_read_xtc(self):
        frame = Frame(
            self.data_dir.joinpath("water.gro"), self.data_dir.joinpath("water.xtc")
        )

        self.check_reference_trajectory(frame)

    def test_frame_accept_path(self):
        """Test that Frame accepts :class:`pathlib.Path` arguments."""
        _ = Frame(
            self.data_dir.joinpath("water.gro"), self.data_dir.joinpath("water.xtc")
        )

    def test_frame_write_xtc(self):
        """Test that :class:`Frame` can save a trajectory file."""
        frame = Frame(
            self.data_dir.joinpath("water.gro"), self.data_dir.joinpath("water.xtc")
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)

            frame.save(tmp_path.joinpath("water.xtc"))
            self.assertTrue(
                util.compare_trajectories(
                    self.data_dir.joinpath("water.xtc"),
                    tmp_path.joinpath("water.xtc"),
                    topology_file=self.data_dir.joinpath("water.gro"),
                )
            )

    def test_raise_nonmatching_system_mdtraj(self):
        with self.assertRaises(NonMatchingSystemError):
            _ = Frame(
                self.data_dir.joinpath("water.gro"), self.data_dir.joinpath("sugar.xtc")
            )


if __name__ == "__main__":
    unittest.main()
