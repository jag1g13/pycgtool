#!/usr/bin/env python3

import argparse
import cProfile
import logging
import pathlib
import sys
import typing

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
    bonds = BondSet(config.bnd, config)

    logger.info("Bond measurements will be made")
    bonds.apply(frame)

    if config.map and config.trajectory:
        # Only perform Boltzmann Inversion if we have a mapping and a trajectory.
        # Otherwise we get infinite force constants.
        logger.info("Beginning Boltzmann Inversion")
        bonds.boltzmann_invert()

        if config.output_forcefield:
            logger.info("Creating GROMACS forcefield directory")
            out_dir = pathlib.Path(config.out_dir)
            forcefield = ForceField(config.output_name, dir_path=out_dir)
            forcefield.write(config.output_name, mapping, bonds)
            logger.info("GROMACS forcefield directory created")

        else:
            bonds.write_itp(get_output_filepath('itp', config),
                            mapping=mapping)

    if config.dump_measurements:
        logger.info("Dumping bond measurements to file")
        bonds.dump_values(config.dump_n_values, config.out_dir)


def mapping_loop(frame: Frame, config) -> typing.Tuple[Frame, Mapping]:
    """Perform mapping loop over input trajectory.

    :param frame:
    :param config: Program arguments from argparse
    """
    logger.info("Mapping will be performed")
    mapping = Mapping(config.map, config, itp_filename=config.itp)

    cg_frame = mapping.apply(frame)
    cg_frame.save(get_output_filepath('gro', config), frame_number=0)

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

    if config.map:
        cg_frame, mapping = mapping_loop(frame, config)

    else:
        logger.info("Mapping will not be performed")
        mapping = None
        cg_frame = frame

    if config.output_xtc:
        cg_frame.save(get_output_filepath('xtc', config))

    if config.bnd:
        measure_bonds(cg_frame, mapping, config)


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
    input_files.add_argument('-m', '--map', type=str,
                             help="Mapping file")
    input_files.add_argument('-b', '--bnd', type=str,
                             help="Bonds file")
    input_files.add_argument('-i', '--itp', type=str,
                             help="GROMACS ITP file")
    input_files.add_argument('--begin', type=int, default=0,
                             help="Frame number to begin")
    input_files.add_argument('--end', type=int, default=None,
                             help="Frame number to end")

    # Output files
    output_files = parser.add_argument_group("output files")

    output_files.add_argument('--out-dir', default='.', type=str,
                              help="Directory where output files should be placed")
    output_files.add_argument('--output-xtc', default=False, action='store_true',
                              help="Output a pseudo-CG trajectory")
    output_files.add_argument("--output-name", default="out",
                              help="Base name of output files")
    output_files.add_argument("--output", default="gro",
                              help="Coordinate output format")
    output_files.add_argument("--output-forcefield", default=False, action='store_true',
                              help="Output GROMACS forefield directory?")
    output_files.add_argument("--dump-measurements", default=False, action="store_true",
                              help="Whether to output bond measurements")
    output_files.add_argument("--dump-n-values", type=int, default=10000,
                              help="How many measurements to output")

    # Mapping options
    mapping_options = parser.add_argument_group("mapping options")

    mapping_options.add_argument("--map-only", default=False, action="store_true",
                                 help="Run in mapping-only mode")
    mapping_options.add_argument("--map-center", default="geom",
                                 choices=["geom", "mass", "first"],
                                 help="Mapping method")
    mapping_options.add_argument("--virtual-map-center", default="geom",
                                 choices=["geom", "mass"],
                                 help="Virtual site mapping method")

    # Bond options
    bond_options = parser.add_argument_group("bond options")

    bond_options.add_argument("--constr_threshold", type=float, default=100000,
                              help="Convert stiff bonds to contraints over [value]")
    bond_options.add_argument("--temperature", type=float, default=310,
                              help="Temperature of reference simulation")
    bond_options.add_argument("--default-fc", default=False, action='store_true',
                              help="Use default MARTINI force constants?")
    bond_options.add_argument("--generate-angles", default=True, action='store_false',
                              help="Generate angles from bonds")
    bond_options.add_argument("--generate-dihedrals", default=False, action="store_true",
                              help="Generate dihedrals from bonds")
    bond_options.add_argument("--length-form", default="Harmonic",
                              help="Form of bond potential")
    bond_options.add_argument("--angle-form", default="CosHarmonic",
                              help="Form of angle potential")
    bond_options.add_argument("--dihedral-form", default="Harmonic",
                              help="Form of dihedral potential")

    # Run options
    run_options = parser.add_argument_group("run options")

    run_options.add_argument('--profile', default=False, action='store_true',
                             help="Profile performance")
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
        args.dump_measurements = bool(args.bnd) and not bool(args.map)

    if not args.map_only:
        args.map_only = not bool(args.bnd)

    if not args.map and not args.bnd:
        raise ArgumentValidationError("One or both of -m and -b is required.")

    return args


def main():
    args = parse_arguments(sys.argv[1:])

    print("Using GRO: {0}".format(args.topology))
    print("Using XTC: {0}".format(args.trajectory))

    if args.profile:
        with cProfile.Profile() as profiler:
            full_run(args)

        profiler.dump_stats('gprof.out')

    else:
        full_run(args)


if __name__ == "__main__":
    main()
