#!/usr/bin/env python3

import argparse
import collections

from pycgtool.frame import Frame
from pycgtool.mapping import Mapping
from pycgtool.bondset import BondSet
from pycgtool.forcefield import ForceField
from pycgtool.util import Progress, Options


def main(args, config):
    """
    Main function of the program PyCGTOOL.

    Performs the complete AA->CG mapping and outputs a GROMACS forcefield directory.

    :param args: Arguments from argparse
    :param config: Configuration dictionary
    """
    frame = Frame(gro=args.gro, xtc=args.xtc)

    if args.bnd:
        bonds = BondSet(args.bnd, config)

    if args.map:
        mapping = Mapping(args.map, config)

    # Main loop - perform mapping and measurement on every frame in XTC
    numframes = frame.numframes if args.frames == -1 else args.frames
    for _ in Progress(numframes, postwhile=frame.next_frame):
        if args.map:
            cgframe = mapping.apply(frame, exclude={"SOL"})
        else:
            cgframe = frame

        if args.bnd:
            bonds.apply(cgframe)

    if args.map:
        cgframe.output("out.gro", format=config.output)

    if args.bnd:
        bonds.boltzmann_invert()
        bonds.write_itp("out.itp", mapping=mapping)
        if config.output_forcefield:
            ff = ForceField("fftest.ff")
            ff.write_rtp("test.rtp", mapping, bonds)

        if config.dump_measurements:
            bonds.dump_values()


def map_only(args, config):
    """
    Perform AA->CG mapping and output coordinate file.

    :param args: Program arguments
    """
    frame = Frame(gro=args.gro)
    mapping = Mapping(args.map, config)
    cgframe = mapping.apply(frame, exclude={"SOL"})
    cgframe.output("out.gro", format=config.output)


def interactive(config):
    """
    Read in options in interactive terminal mode.  Take defaults as options argument.

    :param default_config: Default configuration
    :return: Dictionary of configuration options
    """
    # TODO add interactive terminal mode to select options eg output format instead of CGTOOL config file
    while True:
        line = input(">>").lower()
        if not line:
            break

        if line == "list":
            for key, value in config:
                print("{0:20s} {1}".format(key, value))
            continue

        try:
            key, value = line.split()
            config.set(key, value)
        except KeyError:
            print("Invalid option '{0}' - use 'list' to view options".format(key))
        except ValueError:
            print("Invalid value '{0}' for option '{1}' - use 'list' to view options".format(value, key))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform coarse-grain mapping of atomistic trajectory")
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-g', '--gro', type=str, required=True, help="GROMACS GRO file")
    required.add_argument('-m', '--map', type=str, required=True, help="Mapping file")
    required.add_argument('-x', '--xtc', type=str, help="GROMACS XTC file")
    required.add_argument('-b', '--bnd', type=str, help="Bonds file")

    parser.add_argument('-i', '--interactive', default=False, action='store_true')
    parser.add_argument('-f', '--frames', type=int, default=-1, help="Number of frames to read")

    args = parser.parse_args()
    print("Using GRO: {0}".format(args.gro))
    print("Using XTC: {0}".format(args.xtc))

    config = Options([("output", "gro"),
                      ("map_only", False),
                      ("map_center", "geom"),
                      ("constr_threshold", 100000),
                      ("dump_measurements", False),
                      ("output_forcefield", False)])

    if args.interactive:
        interactive(config)

    if config.map_only or (args.bnd is None and args.xtc is None):
        map_only(args, config)
    else:
        main(args, config)
