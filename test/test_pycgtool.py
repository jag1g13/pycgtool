"""Tests of the full application workflows."""

import itertools
import pathlib
import subprocess
import tempfile
import typing
import unittest

from pycgtool import util
import pycgtool.__main__ as main


def get_args(name, out_dir, extra: typing.Optional[typing.Mapping] = None):
    data_dir = pathlib.Path('test/data')

    args = {
        '-g': data_dir.joinpath(f'{name}.gro'),
        '-x': data_dir.joinpath(f'{name}.xtc'),
        '-m': data_dir.joinpath(f'{name}.map'),
        '-b': data_dir.joinpath(f'{name}.bnd'),
        '--out-dir': out_dir,
    }

    parsed_args = main.parse_arguments(
        itertools.chain(*[[key, str(value)] for key, value in args.items()]))

    if extra is not None:
        for key, value in extra.items():
            setattr(parsed_args, key, value)

    return parsed_args


class PycgtoolTest(unittest.TestCase):
    base_dir = pathlib.Path(__file__).absolute().parent
    data_dir = base_dir.joinpath('data')

    def test_run_help(self):
        self.assertEqual(
            0,
            subprocess.check_call(["python", "-m", "pycgtool", "-h"],
                                  stdout=subprocess.PIPE))

    def test_parse_arguments(self):
        args = main.parse_arguments([
            '-g',
            'GRO',
            '-m',
            'MAP',
            '--begin',
            '1000',
        ])

        self.assertEqual('GRO', args.gro)
        self.assertEqual('MAP', args.map)
        self.assertEqual(1000, args.begin)

    def test_map_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            args = get_args('sugar', tmp_path, extra={
                'output_xtc': True,
                'bnd': None,
            })

            main.full_run(args)

            self.assertTrue(
                util.compare_trajectories(
                    self.data_dir.joinpath('sugar_out.xtc'),
                    tmp_path.joinpath('out.xtc'),
                    topology_file=self.data_dir.joinpath('sugar_out.gro')))

    def test_full(self):
        base_dir = pathlib.Path(__file__).absolute().parent
        data_dir = base_dir.joinpath('data')

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            args = get_args('sugar', tmp_path)

            main.full_run(args)

            self.assertTrue(
                util.cmp_file_whitespace_float(
                    tmp_path.joinpath("out.itp"),
                    data_dir.joinpath("sugar_out.itp"),
                    rtol=0.001,
                    verbose=True))

            self.assertTrue(
                util.compare_trajectories(
                    data_dir.joinpath('sugar_out.gro'),
                    tmp_path.joinpath('out.gro'),
                    topology_file=data_dir.joinpath('sugar_out.gro'),
                    rtol=0.001))

    def test_full_no_traj(self):
        base_dir = pathlib.Path(__file__).absolute().parent
        data_dir = base_dir.joinpath('data')

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            args = get_args('sugar', tmp_path, extra={
                'xtc': None,
            })

            main.full_run(args)

            self.assertFalse(tmp_path.joinpath('out.itp').exists())

            self.assertTrue(
                util.compare_trajectories(
                    data_dir.joinpath('sugar_out.gro'),
                    tmp_path.joinpath('out.gro'),
                    topology_file=data_dir.joinpath('sugar_out.gro'),
                    rtol=0.001))

    # TODO more tests


if __name__ == '__main__':
    unittest.main()
