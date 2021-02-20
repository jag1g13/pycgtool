#!/usr/bin/env python3

import argparse
import cProfile
import logging
import pathlib
import sys
import textwrap
import typing

from mdplus.multiscale import GLIMPS
import numpy as np
from rich.logging import RichHandler

from .frame import Frame
from .mapping import Mapping
from .bondset import BondSet
from .forcefield import ForceField

logger = logging.getLogger(__name__)


class ArgumentValidationError(ValueError):
    """Exception raised for invalid combinations of command line arguments."""


class PyCGTOOL:
    def __init__(self, config):
        self.config = config

    def get_output_filepath(self, ext: str) -> pathlib.Path:
        """Get file path for an output file by extension.

        :param ext:
        """
        out_dir = pathlib.Path(self.config.out_dir)
        return out_dir.joinpath(self.config.output_name + '.' + ext)


    def measure_bonds(self, frame: Frame, mapping: typing.Optional[Mapping]) -> None:
        """Measure bonds at the end of a run.

        :param frame:
        :param mapping:
        """
        bonds = BondSet(self.config.bondset, self.config)
        bonds.apply(frame)

        if self.config.mapping and self.config.trajectory:
            # Only perform Boltzmann Inversion if we have a mapping and a trajectory.
            # Otherwise we get infinite force constants.
            logger.info('Starting Boltzmann Inversion')
            bonds.boltzmann_invert()
            logger.info('Finished Boltzmann Inversion')

            if self.config.output_forcefield:
                logger.info("Writing GROMACS forcefield directory")
                out_dir = pathlib.Path(self.config.out_dir)
                forcefield = ForceField(self.config.output_name, dir_path=out_dir)
                forcefield.write(self.config.output_name, mapping, bonds)
                logger.info("Finished writing GROMACS forcefield directory")

            else:
                bonds.write_itp(self.get_output_filepath('itp'),
                                mapping=mapping)

        if self.config.dump_measurements:
            logger.info('Writing bond measurements to file')
            bonds.dump_values(self.config.dump_n_values, self.config.out_dir)
            logger.info('Finished writing bond measurements to file')


    def mapping_loop(self, frame: Frame) -> typing.Tuple[Frame, Mapping]:
        """Perform mapping loop over input trajectory.

        :param frame:
        """
        logger.info('Starting AA->CG mapping')
        mapping = Mapping(self.config.mapping, self.config, itp_filename=self.config.itp)

        cg_frame = mapping.apply(frame)
        cg_frame.save(self.get_output_filepath('gro'), frame_number=0)
        logging.info('Finished AA->CG mapping')

        self.train_backmapper(frame, cg_frame)

        return cg_frame, mapping


    def get_coords(self, frame: Frame, resname: str) -> np.ndarray:
        return np.concatenate([
            frame._trajectory.atom_slice([atom.index for atom in residue.atoms]).xyz
            for residue in frame._trajectory.topology.residues
            if residue.name == resname
        ])


    def train_backmapper(self, aa_frame: Frame, cg_frame: Frame):
        # resname = 'POPC'
        # aa_coords = get_coords(aa_frame, resname)
        # cg_coords = get_coords(cg_frame, resname)

        cg_subset_traj = cg_frame._trajectory.atom_slice(cg_frame._trajectory.topology.select('resid 1'))
        aa_subset_traj = aa_frame._trajectory.atom_slice(aa_frame._trajectory.topology.select('resid 1'))

        cg_subset_traj.save('cg_test.gro')
        aa_subset_traj.save('aa_test.gro')

        logger.info('Training backmapper')
        backmapper = GLIMPS()
        backmapper.fit(cg_subset_traj.xyz, aa_subset_traj.xyz)
        logger.info('Finished training backmapper')

        logger.info('Testing backmapper')
        backmapped = backmapper.transform(cg_subset_traj.xyz)
        aa_subset_traj.xyz = backmapped
        aa_subset_traj.save('backmapped.gro')
        logger.info('Finished testing backmapper')

        backmapper.save(self.get_output_filepath('backmapper.pkl'))


    def full_run(self):
        """Main function of the program PyCGTOOL.

        Performs the complete AA->CG mapping and outputs a files dependent on given input.
        """
        frame = Frame(
            topology_file=self.config.topology,
            trajectory_file=self.config.trajectory,  # May be None
            frame_start=self.config.begin,
            frame_end=self.config.end)
        frame._trajectory.make_molecules_whole(inplace=True)

        if self.config.mapping:
            cg_frame, mapping = self.mapping_loop(frame)

        else:
            logger.info('Skipping AA->CG mapping')
            mapping = None
            cg_frame = frame

        if self.config.output_xtc:
            cg_frame.save(self.get_output_filepath('xtc'))

        if self.config.bondset:
            self.measure_bonds(cg_frame, mapping)


class BooleanAction(argparse.Action):
    """Set up a boolean argparse argument with matching `--no` argument.

    Based on https://thisdataguy.com/2017/07/03/no-options-with-argparse-and-python/
    """
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        super().__init__(option_strings, dest, nargs=0, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, not option_string.startswith('--no-'))


def parse_arguments(arg_list):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=
        "Generate coarse-grained molecular dynamics models from atomistic trajectories."
    )

    # yapf: disable
    # Input files
    input_files = parser.add_argument_group("input files")

    input_files.add_argument('topology', type=str,
                             help="AA simulation topology - e.g. PDB, GRO, etc.")
    input_files.add_argument('trajectory', type=str, nargs='?',
                             help="AA simulation trajectory - e.g. XTC, DCD, etc.")
    input_files.add_argument('-m', '--mapping', type=str,
                             help="Mapping file")
    input_files.add_argument('-b', '--bondset', type=str,
                             help="Bonds file")
    input_files.add_argument('-i', '--itp', type=str,
                             help="GROMACS ITP file")
    input_files.add_argument('--begin', type=int, default=0,
                             help="Trajectory frame number to begin")
    input_files.add_argument('--end', type=int, default=None,
                             help="Trajectory frame number to end")

    # Output files
    output_files = parser.add_argument_group("output files")

    output_files.add_argument('--out-dir', default='.', type=str,
                              help="Directory where output files should be placed")
    output_files.add_argument("--output-name", default="out",
                              help="Base name of output files")
    output_files.add_argument('--output-xtc', '--no-output-xtc', default=False, action=BooleanAction,
                              help="Output a pseudo-CG trajectory?")
    output_files.add_argument("--output", default="gro",
                              help="Coordinate output format")
    output_files.add_argument("--output-forcefield", '--no-output-forcefield', default=False, action=BooleanAction,
                              help="Output GROMACS forefield directory?")
    output_files.add_argument("--dump-measurements", '--no-dump-measurements', default=False, action=BooleanAction,
                              help="Output sample of bond measurements?")
    output_files.add_argument("--dump-n-values", type=int, default=10000,
                              help="Size of sample of measurements to output")

    # Mapping options
    mapping_options = parser.add_argument_group("mapping options")

    mapping_options.add_argument("--map-only", '--no-map-only', default=False, action=BooleanAction,
                                 help="Run in mapping-only mode?")
    mapping_options.add_argument("--map-center", default="geom",
                                 choices=["geom", "mass", "first"],
                                 help="Mapping method")
    mapping_options.add_argument("--virtual-map-center", default="geom",
                                 choices=["geom", "mass"],
                                 help="Virtual site mapping method")

    # Bond options
    bond_options = parser.add_argument_group("bond options")

    bond_options.add_argument("--constr_threshold", type=float, default=100000,
                              help="Convert bonds with force constants over [value] to constraints")
    bond_options.add_argument("--temperature", type=float, default=310,
                              help="Temperature of reference simulation")
    bond_options.add_argument("--default-fc", '--no-default-fc', default=False, action=BooleanAction,
                              help="Use default MARTINI force constants?")
    bond_options.add_argument("--generate-angles", '--no-generate-angles', default=True, action=BooleanAction,
                              help="Generate angles from bonds?")
    bond_options.add_argument("--generate-dihedrals", '--no-generate-dihedrals', default=False, action=BooleanAction,
                              help="Generate dihedrals from bonds?")
    bond_options.add_argument("--length-form", default="Harmonic",
                              help="Form of bond potential")
    bond_options.add_argument("--angle-form", default="CosHarmonic",
                              help="Form of angle potential")
    bond_options.add_argument("--dihedral-form", default="Harmonic",
                              help="Form of dihedral potential")

    # Run options
    run_options = parser.add_argument_group("run options")

    run_options.add_argument('--profile', '--no-profile', default=False, action=BooleanAction,
                             help="Profile performance?")
    run_options.add_argument('--log-level', default='INFO',
                             choices=('DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'),
                             help="Which log messages should be shown?")
    # yapf: enable

    args = parser.parse_args(arg_list)

    try:
        args = validate_arguments(args)

    except ArgumentValidationError as exc:
        parser.error(exc)

    return args


def validate_arguments(args):
    """Check that arguments are not contradictory and modify where necessary.

    :param args: Parsed arguments from ArgumentParser
    """
    if not args.dump_measurements:
        args.dump_measurements = bool(args.bondset) and not bool(args.mapping)

    if not args.map_only:
        args.map_only = not bool(args.bondset)

    if not args.mapping and not args.bondset:
        raise ArgumentValidationError("One or both of -m and -b is required.")

    return args


def main():
    args = parse_arguments(sys.argv[1:])

    logging.basicConfig(level=args.log_level,
                        format='%(message)s',
                        datefmt='[%X]',
                        handlers=[RichHandler(rich_tracebacks=True)])

    banner = """\
         _____        _____ _____ _______ ____   ____  _      
        |  __ \      / ____/ ____|__   __/ __ \ / __ \| |     
        | |__) |   _| |   | |  __   | | | |  | | |  | | |     
        |  ___/ | | | |   | | |_ |  | | | |  | | |  | | |     
        | |   | |_| | |___| |__| |  | | | |__| | |__| | |____ 
        |_|    \__, |\_____\_____|  |_|  \____/ \____/|______|
                __/ |                                         
               |___/                                          
    """  # noqa

    logger.info('[bold blue]%s[/]',
                textwrap.dedent(banner),
                extra={'markup': True})

    logger.info(30 * '-')
    logger.info('Topology:\t%s', args.topology)
    logger.info('Trajectory:\t%s', args.trajectory)
    logger.info('Mapping:\t%s', args.mapping)
    logger.info('Bondset:\t%s', args.bondset)
    logger.info(30 * '-')

    try:
        pycgtool = PyCGTOOL(args)

        if args.profile:
            with cProfile.Profile() as profiler:
                pycgtool.full_run()

            profiler.dump_stats('gprof.out')

        else:
            pycgtool.full_run()

        logger.info('Finished processing - goodbye!')

    except Exception as exc:
        logger.exception(exc)


if __name__ == "__main__":
    main()
