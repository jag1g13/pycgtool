#!/usr/bin/env python3

import argparse
import logging
import pathlib
import sys

from .frame import Frame
from .mapping import Mapping
from .bondset import BondSet
from .forcefield import ForceField
from .functionalforms import FunctionalForms
from .interface import Options, Progress

logger = logging.getLogger(__name__)


def full_run(args, config):
    """
    Main function of the program PyCGTOOL.

    Performs the complete AA->CG mapping and outputs a files dependent on given input.

    :param args: Arguments from argparse
    :param config: Configuration dictionary
    """
    frame = Frame(gro=args.gro,
                  xtc=args.xtc,
                  itp=args.itp,
                  frame_start=args.begin)

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

        mapping = Mapping(args.map, config, itp=args.itp)
        cgframe = mapping.apply(frame)
        cgframe.output(out_gro_file, format=config.output)

    else:
        logger.info("Mapping will not be performed")
        cgframe = frame

    # Only measure bonds from GRO frame if no XTC is provided
    # Allows the user to get a topology from a single snapshot
    if args.bnd and args.xtc is None:
        bonds.apply(cgframe)

    # Main loop - perform mapping and measurement on every frame in XTC
    def main_loop():
        nonlocal cgframe
        if not frame.next_frame():
            return False

        if args.map:
            cgframe = mapping.apply(frame, cgframe=cgframe)

            if config.output_xtc:
                cgframe.write_xtc(out_xtc_file)

        else:
            cgframe = frame

        if args.bnd:
            bonds.apply(cgframe)

        return True

    numframes = frame.numframes - args.begin if args.end == -1 else args.end - args.begin
    logger.info("Beginning analysis of {0} frames".format(numframes))
    Progress(numframes, dowhile=main_loop, quiet=args.quiet).run()

    if args.bnd:
        if args.map:
            logger.info("Beginning Boltzmann inversion")
            bonds.boltzmann_invert(progress=(not args.quiet))

            if config.output_forcefield:
                logger.info("Creating GROMACS forcefield directory")
                ff = ForceField(config.output_name, dir_path=out_dir)
                ff.write(config.output_name, mapping, bonds)
                logger.info("GROMACS forcefield directory created")

            else:
                bonds.write_itp(out_itp_file, mapping=mapping)

        if config.dump_measurements:
            logger.info("Dumping bond measurements to file")
            bonds.dump_values(config.dump_n_values)


def map_only(args, config):
    """
    Perform AA->CG mapping and output coordinate file.

    :param args: Program arguments
    :param config: Object containing run options
    """
    frame = Frame(gro=args.gro, xtc=args.xtc)
    mapping = Mapping(args.map, config)

    out_dir = pathlib.Path(args.out_dir)
    out_gro_file = out_dir.joinpath(config.output_name + ".gro")
    out_xtc_file = out_dir.joinpath(config.output_name + ".xtc")

    cgframe = mapping.apply(frame)
    cgframe.output(out_gro_file, format=config.output)

    if args.xtc and (config.output_xtc or args.outputxtc):
        # Main loop - perform mapping and measurement on every frame in XTC
        def main_loop():
            nonlocal cgframe
            if not frame.next_frame():
                return False

            cgframe = mapping.apply(frame, cgframe=cgframe)
            cgframe.write_xtc(out_xtc_file)
            return True

        numframes = frame.numframes - args.begin if args.end == -1 else args.end - args.begin
        logger.info("Beginning analysis of {0} frames".format(numframes))

        # Run main loop with progress bar - ignore returned value
        _ = Progress(numframes, dowhile=main_loop, quiet=args.quiet).run()


def main():
    parser = argparse.ArgumentParser(
        description="Perform coarse-grain mapping of atomistic trajectory")
    input_files = parser.add_argument_group("Input files")

    # yapf: disable
    input_files.add_argument('-g', '--gro', type=str, required=True, help="GROMACS GRO file")
    input_files.add_argument('-m', '--map', type=str, help="Mapping file")
    input_files.add_argument('-x', '--xtc', type=str, help="GROMACS XTC file")
    input_files.add_argument('-b', '--bnd', type=str, help="Bonds file")
    input_files.add_argument('-i', '--itp', type=str, help="GROMACS ITP file")
    # yapf: enable

    parser.add_argument('--advanced',
                        default=False,
                        action='store_true',
                        help="Show advanced options menu")
    parser.add_argument('--out-dir',
                        default='.',
                        type=str,
                        help="Directory where output files should be placed")
    parser.add_argument('--outputxtc',
                        default=False,
                        action='store_true',
                        help="Output a pseudo-CG trajectory")
    parser.add_argument('--quiet',
                        default=False,
                        action='store_true',
                        help="Hide progress bars")
    input_files.add_argument('--begin',
                             type=int,
                             default=0,
                             help="Frame number to begin")
    input_files.add_argument('--end',
                             type=int,
                             default=-1,
                             help="Frame number to end")

    # Populate functional forms dictionary - ignore returned value
    _ = FunctionalForms()

    args = parser.parse_args()
    config = Options(
        [("output_name", "out"), ("output", "gro"),
         ("output_xtc", args.outputxtc), ("map_only", not bool(args.bnd)),
         ("map_center", "geom"), ("constr_threshold", 100000),
         ("dump_measurements", bool(args.bnd) and not bool(args.map)),
         ("dump_n_values", 10000), ("output_forcefield", False),
         ("temperature", 310), ("default_fc", False),
         ("generate_angles", True), ("generate_dihedrals", False),
         ("length_form", "harmonic"), ("angle_form", "cosharmonic"),
         ("dihedral_form", "harmonic")], args)

    if not args.map and not args.bnd:
        parser.error("One or both of -m and -b is required.")

    if args.advanced:
        try:
            config.interactive()

        except KeyboardInterrupt:
            sys.exit(0)

    else:
        print("Using GRO: {0}".format(args.gro))
        print("Using XTC: {0}".format(args.xtc))

    if config.map_only:
        map_only(args, config)

    else:
        full_run(args, config)


if __name__ == "__main__":
    main()
