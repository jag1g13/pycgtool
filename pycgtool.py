#!/usr/bin/env python3

import argparse
import sys
import logging

from pycgtool.frame import Frame
from pycgtool.mapping import Mapping
from pycgtool.bondset import BondSet
from pycgtool.forcefield import ForceField
from pycgtool.interface import Progress, Options

logger = logging.getLogger(__name__)


def main(args, config):
    """
    Main function of the program PyCGTOOL.

    Performs the complete AA->CG mapping and outputs a files dependent on given input.

    :param args: Arguments from argparse
    :param config: Configuration dictionary
    """
    frame = Frame(gro=args.gro, xtc=args.xtc, itp=args.itp, frame_start=args.begin, xtc_reader="mdtraj")

    if args.bnd:
        bonds = BondSet(args.bnd, config)
        logger.info("Bond measurements will be made")
    else:
        logger.info("Bond measurements will not be made")

    if args.map:
        mapping = Mapping(args.map, config, itp=args.itp)
        cgframe = mapping.apply(frame, exclude={"SOL"})
        cgframe.output(config.output_name + ".gro", format=config.output)
        logger.info("Mapping will be performed")
    else:
        logger.info("Mapping will not be performed")

    # Main loop - perform mapping and measurement on every frame in XTC
    def main_loop():
        nonlocal cgframe
        frame.next_frame()
        if args.map:
            cgframe = mapping.apply(frame, cgframe=cgframe, exclude={"SOL"})
            if config.output_xtc:
                cgframe.write_xtc(config.output_name + ".xtc")
        else:
            cgframe = frame
        if args.bnd:
            bonds.apply(cgframe)

    numframes = frame.numframes - args.begin if args.end == -1 else args.end - args.begin
    logger.info("Beginning analysis of {0} frames".format(numframes))
    Progress(numframes, postwhile=main_loop).run()

    if args.bnd:
        if args.map:
            logger.info("Beginning Boltzmann inversion")
            bonds.boltzmann_invert()
            if config.output_forcefield:
                logger.info("Creating GROMACS forcefield directory")
                ff = ForceField("ff" + config.output_name + ".ff")
                ff.write_rtp(config.output_name + ".rtp", mapping, bonds)
                logger.info("GROMACS forcefield directory created")
            else:
                bonds.write_itp(config.output_name + ".itp", mapping=mapping)

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
    cgframe = mapping.apply(frame, exclude={"SOL"})
    cgframe.output(config.output_name + ".gro", format=config.output)

    if args.xtc and (config.output_xtc or args.outputxtc):
        # Main loop - perform mapping and measurement on every frame in XTC
        def main_loop():
            nonlocal cgframe
            frame.next_frame()
            cgframe = mapping.apply(frame, cgframe=cgframe, exclude={"SOL"})
            cgframe.write_xtc(config.output_name + ".xtc")

        numframes = frame.numframes - args.begin if args.end == -1 else args.end - args.begin
        logger.info("Beginning analysis of {0} frames".format(numframes))
        Progress(numframes, postwhile=main_loop).run()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform coarse-grain mapping of atomistic trajectory")
    input_files = parser.add_argument_group("Input files")
    input_files.add_argument('-g', '--gro', type=str, required=True, help="GROMACS GRO file")
    input_files.add_argument('-m', '--map', type=str, help="Mapping file")
    input_files.add_argument('-x', '--xtc', type=str, help="GROMACS XTC file")
    input_files.add_argument('-b', '--bnd', type=str, help="Bonds file")
    input_files.add_argument('-i', '--itp', type=str, help="GROMACS ITP file")

    parser.add_argument('--interactive', default=False, action='store_true')
    parser.add_argument('--outputxtc', default=False, action='store_true')
    input_files.add_argument('--begin', type=int, default=0, help="Frame number to begin")
    input_files.add_argument('--end', type=int, default=-1, help="Frame number to end")

    args = parser.parse_args()
    config = Options([("output_name", "out"),
                      ("output", "gro"),
                      ("output_xtc", args.outputxtc),
                      ("map_only", not bool(args.bnd)),
                      ("map_center", "geom"),
                      ("constr_threshold", 100000),
                      ("dump_measurements", bool(args.bnd) and not bool(args.map)),
                      ("dump_n_values", 100000),
                      ("output_forcefield", False),
                      ("temperature", 310),
                      ("angle_default_fc", True),
                      ("generate_angles", True),
                      ("generate_dihedrals", False)],
                     args)

    if not args.map and not args.bnd:
        parser.error("One or both of -m and -b is required.")

    if args.interactive:
        try:
            config.interactive()
        except KeyboardInterrupt:
            sys.exit(0)
    else:
        print("Using GRO: {0}".format(args.gro))
        print("Using XTC: {0}".format(args.xtc))

    if config.map_only:
        logger.info("Starting: mapping only")
        map_only(args, config)
    else:
        logger.info("Starting: parameter measurement")
        main(args, config)

    logger.info("Analysis complete")
