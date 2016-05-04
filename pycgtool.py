#!/usr/bin/env python3

import argparse

from pycgtool.frame import Frame
from pycgtool.mapping import Mapping
from pycgtool.bondset import BondSet
from pycgtool.forcefield import ForceField
from pycgtool.interface import Progress, Options


def main(args, config):
    """
    Main function of the program PyCGTOOL.

    Performs the complete AA->CG mapping and outputs a files dependent on given input.

    :param args: Arguments from argparse
    :param config: Configuration dictionary
    """
    frame = Frame(gro=args.gro, xtc=args.xtc, itp=args.itp, frame_start=args.begin)

    if args.bnd:
        bonds = BondSet(args.bnd, config)

    if args.map:
        mapping = Mapping(args.map, config, itp=args.itp)
        cgframe = mapping.apply(frame, exclude={"SOL"})
        cgframe.output("out.gro", format=config.output)

    # Main loop - perform mapping and measurement on every frame in XTC
    numframes = frame.numframes - args.begin if args.end == -1 else args.end - args.begin
    for _ in Progress(numframes, postwhile=frame.next_frame):
        if args.map:
            cgframe = mapping.apply(frame, exclude={"SOL"})
        else:
            cgframe = frame

        if args.bnd:
            bonds.apply(cgframe)

    if args.bnd:
        if args.map:
            bonds.boltzmann_invert()
            bonds.write_itp("out.itp", mapping=mapping)
            if config.output_forcefield:
                ff = ForceField("fftest.ff")
                ff.write_rtp("test.rtp", mapping, bonds)

        if config.dump_measurements:
            bonds.dump_values(config.dump_n_values)


def map_only(args, config):
    """
    Perform AA->CG mapping and output coordinate file.

    :param args: Program arguments
    :param config: Object containing run options
    """
    frame = Frame(gro=args.gro)
    mapping = Mapping(args.map, config)
    cgframe = mapping.apply(frame, exclude={"SOL"})
    cgframe.output("out.gro", format=config.output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform coarse-grain mapping of atomistic trajectory")
    input_files = parser.add_argument_group("Input files")
    input_files.add_argument('-g', '--gro', type=str, required=True, help="GROMACS GRO file")
    input_files.add_argument('-m', '--map', type=str, help="Mapping file")
    input_files.add_argument('-x', '--xtc', type=str, help="GROMACS XTC file")
    input_files.add_argument('-b', '--bnd', type=str, help="Bonds file")
    input_files.add_argument('-i', '--itp', type=str, help="GROMACS ITP file")

    parser.add_argument('--interactive', default=False, action='store_true')
    # parser.add_argument('-f', '--frames', type=int, default=-1, help="Number of frames to read")
    input_files.add_argument('--begin', type=int, default=0, help="Frame number to begin")
    input_files.add_argument('--end', type=int, default=-1, help="Frame number to end")

    args = parser.parse_args()
    config = Options([("output", "gro"),
                      ("map_only", False),
                      ("map_center", "geom"),
                      ("constr_threshold", 100000),
                      ("dump_measurements", False),
                      ("dump_n_values", 100000),
                      ("output_forcefield", False),
                      ("temperature", 310),
                      ("angle_default_fc", True),
                      ("generate_angles", True),
                      ("generate_dihedrals", False)],
                     args)
    if not args.bnd:
        config.set("map_only", True)

    if args.bnd and not args.map:
        config.set("dump_measurements", True)

    if not args.map and not args.bnd:
        parser.error("One or both of -m and -b is required.")

    if args.interactive:
        config.interactive()
    else:
        print("Using GRO: {0}".format(args.gro))
        print("Using XTC: {0}".format(args.xtc))

    if config.map_only:
        map_only(args, config)
    else:
        main(args, config)
