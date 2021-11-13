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
    base_dir = pathlib.Path(__file__).absolute().parent
    data_dir = base_dir.joinpath('data')

    args = [
        data_dir.joinpath(f'{name}.gro'),
        data_dir.joinpath(f'{name}.xtc'),
    ]

    kwargs = {
        '-m': data_dir.joinpath(f'{name}.map'),
        '-b': data_dir.joinpath(f'{name}.bnd'),
        '--out-dir': out_dir,
    }

    parsed_args = main.parse_arguments(
        itertools.chain(map(str, args),
                        *[[key, str(value)] for key, value in kwargs.items()]))

    if extra is not None:
        for key, value in extra.items():
            setattr(parsed_args, key, value)

    # Re-validate after manual changes
    return main.validate_arguments(parsed_args)


class PycgtoolTest(unittest.TestCase):
    base_dir = pathlib.Path(__file__).absolute().parent
    data_dir = base_dir.joinpath('data')

    def test_pycgtool_help(self):
        subprocess.check_call(["python", "-m", "pycgtool", "-h"])
        subprocess.check_call(["python", "-m", "pycgtool", "--help"])

    def test_pycgtool_version(self):
        subprocess.check_call(['python', '-m', 'pycgtool', '-v'])
        subprocess.check_call(['python', '-m', 'pycgtool', '--version'])

    def test_parse_arguments(self):
        args = main.parse_arguments([
            'TOPOLOGY',
            '-m',
            'MAP',
            '--begin',
            '1000',
        ])

        self.assertEqual('TOPOLOGY', args.topology)
        self.assertEqual('MAP', args.mapping)
        self.assertEqual(1000, args.begin)

    def test_map_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            args = get_args('sugar', tmp_path, extra={
                'output_xtc': True,
                'bnd': None,
            })

            # Equivalent to
            # pycgtool <top> <trj> -m <map> --output-xtc
            main.PyCGTOOL(args)

            self.assertTrue(
                util.compare_trajectories(
                    self.data_dir.joinpath('sugar_out.xtc'),
                    tmp_path.joinpath('out.xtc'),
                    topology_file=self.data_dir.joinpath('sugar_out.gro')))

    def test_full(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = pathlib.Path(tmpdir)
            args = get_args('sugar', tmp_path)

            # Equivalent to
            # pycgtool <top> <trj> -m <map> -b <bnd>
            main.PyCGTOOL(args)

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
                'trajectory': None,
            })

            # Equivalent to
            # pycgtool <top> -m <map> -b <bnd>
            main.PyCGTOOL(args)

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
                'mapping': None,
                'trajectory': None,
            })

            # Equivalent to
            # pycgtool <top> -b <bnd>
            main.PyCGTOOL(args)

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

            # Equivalent to
            # pycgtool <top> <trj> -m <map> -b <bnd> --output-forcefield
            main.PyCGTOOL(args)

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
