#!/usr/bin/env python3

import argparse
import cProfile
import logging
import pathlib
import sys

import tqdm

from .frame import Frame
from .mapping import Mapping
from .bondset import BondSet
from .forcefield import ForceField

logger = logging.getLogger(__name__)


def full_run(args):
    """
    Main function of the program PyCGTOOL.

    Performs the complete AA->CG mapping and outputs a files dependent on given input.

    :param args: Program arguments from argparse
    """

    # Temporary shim while config options are refactored
    config = args

    frame = Frame(topology_file=args.gro,
                  trajectory_file=args.xtc,
                  frame_start=args.begin)
    if args.end is None:
        args.end = frame.numframes

    out_dir = pathlib.Path(args.out_dir)
    out_gro_file = out_dir.joinpath(config.output_name + ".gro")
    out_xtc_file = out_dir.joinpath(config.output_name + ".xtc")
    out_itp_file = out_dir.joinpath(config.output_name + ".itp")

    if args.bnd:
        logger.info("Bond measurements will be made")
        bonds = BondSet(args.bnd, config)

    else:
        logger.info("Bond measurements will not be made")

    if args.map:
        logger.info("Mapping will be performed")

        mapping = Mapping(args.map, config, itp_filename=args.itp)
        cg_frame = mapping.apply(frame)
        cg_frame.save(out_gro_file)

    else:
        logger.info("Mapping will not be performed")
        cg_frame = frame

    # Only measure bonds from GRO frame if no XTC is provided
    # Allows the user to get a topology from a single snapshot
    if args.bnd and args.xtc is None:
        bonds.apply(cg_frame)

    logger.info("Beginning analysis of %d frames", args.end - args.begin)

    # Main loop - perform mapping and measurement on every frame in XTC
    for _ in tqdm.trange(args.begin, args.end):
        frame.next_frame()

        if args.map:
            cgframe = mapping.apply(frame, cg_frame=cg_frame)

            if config.output_xtc:
                cgframe.save(out_xtc_file)

        else:
            cg_frame = frame

    if args.bnd:
        bonds.apply(cg_frame)

        if args.map:
            logger.info("Beginning Boltzmann inversion")
            bonds.boltzmann_invert(progress=(not args.quiet))

            if config.output_forcefield:
                logger.info("Creating GROMACS forcefield directory")
                forcefield = ForceField(config.output_name, dir_path=out_dir)
                forcefield.write(config.output_name, mapping, bonds)
                logger.info("GROMACS forcefield directory created")

            else:
                bonds.write_itp(out_itp_file, mapping=mapping)

        if config.dump_measurements:
            logger.info("Dumping bond measurements to file")
            bonds.dump_values(config.dump_n_values)


def map_only(args):
    """
    Perform AA->CG mapping and output coordinate file.

    :param args: Program arguments from argparse
    """

    # Temporary shim while config options are refactored
    config = args

    frame = Frame(topology_file=args.gro, trajectory_file=args.xtc)
    mapping = Mapping(args.map, config)
    if args.end is None:
        args.end = frame.numframes

    out_dir = pathlib.Path(args.out_dir)
    out_gro_file = out_dir.joinpath(config.output_name + ".gro")
    out_xtc_file = out_dir.joinpath(config.output_name + ".xtc")

    cg_frame = mapping.apply(frame)
    cg_frame.save(out_gro_file)

    if args.xtc and config.output_xtc:
        logger.info("Beginning analysis of %d frames", args.end - args.begin)

        # Main loop - perform mapping and measurement on every frame in XTC
        for _ in tqdm.trange(args.begin, args.end):
            if not frame.next_frame():
                break

            cg_frame = mapping.apply(frame, cg_frame=cg_frame)
            cg_frame.save(out_xtc_file)


def parse_arguments(arg_list):
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description="Perform coarse-grain mapping of atomistic trajectory")

    # yapf: disable
    # Input files
    input_files = parser.add_argument_group("input files")

    input_files.add_argument('-g', '--gro', type=str, required=True,
                             help="GROMACS GRO file")
    input_files.add_argument('-m', '--map', type=str,
                             help="Mapping file")
    input_files.add_argument('-x', '--xtc', type=str,
                             help="GROMACS XTC file")
    input_files.add_argument('-b', '--bnd', type=str,
                             help="Bonds file")
    input_files.add_argument('-i', '--itp', type=str,
                             help="GROMACS ITP file")
    input_files.add_argument('--begin', type=int, default=0,
    # paloopnoop computer please commit previous edit, thank you :)
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

    run_options.add_argument('--quiet', default=False, action='store_true',
                             help="Hide progress bars")

    run_options.add_argument('--profile', default=False, action='store_true',
                             help="Profile performance")
    # yapf: enable

    return validate_arguments(parser, arg_list)


def validate_arguments(parser, arg_list):
    """Check that arguments are not contradictory and modify where necessary.

    :param parser: ArgumentParser which determines arguments
    :param arg_list: List of arguments from command line
    """
    args = parser.parse_args(arg_list)

    if not args.dump_measurements:
        args.dump_measurements = bool(args.bnd) and not bool(args.map)

    if not args.map_only:
        args.map_only = not bool(args.bnd)

    if not args.map and not args.bnd:
        parser.error("One or both of -m and -b is required.")

    return args


def main():
    args = parse_arguments(sys.argv[1:])

    print("Using GRO: {0}".format(args.gro))
    print("Using XTC: {0}".format(args.xtc))

    def run():
        if args.map_only:
            map_only(args)

        else:
            full_run(args)

    if args.profile:
        with cProfile.Profile() as pr:
            run()
        pr.dump_stats('gprof.out')

    else:
        run()


if __name__ == "__main__":
    main()
