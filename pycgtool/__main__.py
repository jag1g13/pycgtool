#!/usr/bin/env python3

import argparse
import cProfile
import logging
import pathlib
import sys
import textwrap
import typing

from rich.logging import RichHandler

from .frame import Frame
from .mapping import Mapping
from .bondset import BondSet
from .forcefield import ForceField

logger = logging.getLogger(__name__)


class ArgumentValidationError(ValueError):
    """Exception raised for invalid combinations of command line arguments."""


def get_output_filepath(ext: str, config) -> pathlib.Path:
    """Get file path for an output file by extension.

    :param ext:
    :param config: Program arguments from argparse
    """
    out_dir = pathlib.Path(config.out_dir)
    return out_dir.joinpath(config.output_name + '.' + ext)


def measure_bonds(frame: Frame, mapping: typing.Optional[Mapping],
                  config) -> None:
    """Measure bonds at the end of a run.

    :param frame:
    :param mapping:
    :param config: Program arguments from argparse
    """
    bonds = BondSet(config.bondset, config)
    bonds.apply(frame)

    if config.mapping and config.trajectory:
        # Only perform Boltzmann Inversion if we have a mapping and a trajectory.
        # Otherwise we get infinite force constants.
        logger.info('Starting Boltzmann Inversion')
        bonds.boltzmann_invert()
        logger.info('Finished Boltzmann Inversion')

        if config.output_forcefield:
            logger.info("Writing GROMACS forcefield directory")
            out_dir = pathlib.Path(config.out_dir)
            forcefield = ForceField(config.output_name, dir_path=out_dir)
            forcefield.write(config.output_name, mapping, bonds)
            logger.info("Finished writing GROMACS forcefield directory")

        else:
            bonds.write_itp(get_output_filepath('itp', config),
                            mapping=mapping)

    if config.dump_measurements:
        logger.info('Writing bond measurements to file')
        bonds.dump_values(config.dump_n_values, config.out_dir)
        logger.info('Finished writing bond measurements to file')


def mapping_loop(frame: Frame, config) -> typing.Tuple[Frame, Mapping]:
    """Perform mapping loop over input trajectory.

    :param frame:
    :param config: Program arguments from argparse
    """
    logger.info('Starting AA->CG mapping')
    mapping = Mapping(config.mapping, config, itp_filename=config.itp)

    cg_frame = mapping.apply(frame)
    cg_frame.save(get_output_filepath('gro', config), frame_number=0)
    logging.info('Finished AA->CG mapping')

    return cg_frame, mapping


def full_run(config):
    """Main function of the program PyCGTOOL.

    Performs the complete AA->CG mapping and outputs a files dependent on given input.

    :param config: Program arguments from argparse
    """
    frame = Frame(
        topology_file=config.topology,
        trajectory_file=config.trajectory,  # May be None
        frame_start=config.begin,
        frame_end=config.end)

    if config.mapping:
        cg_frame, mapping = mapping_loop(frame, config)

    else:
        logger.info('Skipping AA->CG mapping')
        mapping = None
        cg_frame = frame

    if config.output_xtc:
        cg_frame.save(get_output_filepath('xtc', config))

    if config.bondset:
        measure_bonds(cg_frame, mapping, config)


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
                        handlers=[RichHandler()])

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
        if args.profile:
            with cProfile.Profile() as profiler:
                full_run(args)

            profiler.dump_stats('gprof.out')

        else:
            full_run(args)

        logger.info('Finished processing - goodbye!')

    except Exception as exc:
        logger.error(exc)


if __name__ == "__main__":
    main()
