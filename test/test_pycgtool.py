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

    # Re-validate after manual changes
    return main.validate_arguments(parsed_args)


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
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            args = get_args('sugar', tmp_path)

            main.full_run(args)

            self.assertTrue(
                util.cmp_file_whitespace_float(
                    tmp_path.joinpath("out.itp"),
                    self.data_dir.joinpath("sugar_out.itp"),
                    rtol=0.001,
                    verbose=True))

            self.assertTrue(
                util.compare_trajectories(
                    self.data_dir.joinpath('sugar_out.gro'),
                    tmp_path.joinpath('out.gro'),
                    rtol=0.001))

    def test_full_no_traj(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            args = get_args('sugar', tmp_path, extra={
                'xtc': None,
            })

            main.full_run(args)

            self.assertFalse(tmp_path.joinpath('out.itp').exists())

            self.assertTrue(
                util.compare_trajectories(
                    self.data_dir.joinpath('sugar_out.gro'),
                    tmp_path.joinpath('out.gro'),
                    rtol=0.001))

    def test_measure_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            args = get_args('sugar', tmp_path, extra={
                'map': None,
                'xtc': None,
            })

            main.full_run(args)

            # Does not produce itp file
            self.assertFalse(tmp_path.joinpath('out.itp').exists())

            # Check bond dump files against reference
            for bond_type in ['length', 'angle', 'dihedral']:
                out_file = tmp_path.joinpath(f'ALLA_{bond_type}.dat')
                self.assertTrue(out_file.exists())

                self.assertTrue(util.cmp_file_whitespace_float(
                    self.data_dir.joinpath(f'ALLA_{bond_type}_one.dat'),
                    out_file
                ))

    def test_forcefield(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            args = get_args('sugar',
                            tmp_path,
                            extra={
                                'output_forcefield': True,
                            })

            main.full_run(args)

            # Does not produce itp file
            self.assertFalse(tmp_path.joinpath('out.itp').exists())

            out_ff_dir = tmp_path.joinpath('ffout.ff')
            self.assertTrue(out_ff_dir.is_dir())

            # Compare all files in ffout.ff to reference versions
            for out_file in out_ff_dir.iterdir():
                ref_file = self.data_dir.joinpath(out_file)

                self.assertTrue(
                    util.cmp_file_whitespace_float(ref_file, out_file))

    # TODO more tests
    # TODO test wrong args.end


if __name__ == '__main__':
    unittest.main()
