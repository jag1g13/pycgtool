import os
import pathlib
import subprocess
import tempfile
import unittest

from pycgtool import util
from pycgtool.__main__ import map_only


class Args:
    itp = None
    begin = 0
    end = -1
    quiet = True

    map_center = "geom"
    virtual_map_center = "geom"
    output_xtc = True
    output_name = "out"
    output = "gro"

    def __init__(self, name, use_map_file=True, use_bnd_file=True):
        self.gro = os.path.join("test/data", name + ".gro")
        self.xtc = os.path.join("test/data", name + ".xtc")
        self.map = os.path.join("test/data", name + ".map") if use_map_file else None
        self.bnd = os.path.join("test/data", name + ".bnd") if use_bnd_file else None


class PycgtoolTest(unittest.TestCase):
    def test_run_help(self):
        self.assertEqual(
            0,
            subprocess.check_call(["python", "-m", "pycgtool", "-h"],
                                  stdout=subprocess.PIPE))

    def test_map_only(self):
        map_args = Args("sugar")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            map_args.out_dir = tmp_path  # pylint: disable=attribute-defined-outside-init

            map_only(map_args)

            self.assertTrue(util.compare_trajectories(
                'test/data/sugar_out.xtc',
                tmp_path.joinpath('out.xtc'),
                topology_file='test/data/sugar_out.gro'
            ))

    def test_full(self):
        base_dir = pathlib.Path(__file__).absolute().parent
        data_dir = base_dir.joinpath('data')

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            return_code = subprocess.check_call([
                "python", "-m", "pycgtool",
                "-g", data_dir.joinpath("sugar.gro"),
                "-x", data_dir.joinpath("sugar.xtc"),
                "-m", data_dir.joinpath("sugar_only.map"),
                "-b", data_dir.joinpath("sugar.bnd"),
                "--out-dir", tmpdir,
            ])  # yapf: disable

            self.assertEqual(0, return_code)

            self.assertTrue(
                util.cmp_file_whitespace_float(tmp_path.joinpath("out.itp"),
                                          data_dir.joinpath("sugar_out.itp"),
                                          rtol=0.001,
                                          verbose=True))
            self.assertTrue(util.compare_trajectories(
                'test/data/sugar_out.gro',
                tmp_path.joinpath('out.gro'),
                topology_file='test/data/sugar_out.gro',
                rtol=0.001
            ))

    # TODO more tests


if __name__ == '__main__':
    unittest.main()
