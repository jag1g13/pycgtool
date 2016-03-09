#!/usr/bin/env python3

import argparse

from pycgtool.frame import Frame
from pycgtool.mapping import Mapping
from pycgtool.bondset import BondSet
from pycgtool.forcefield import ForceField


def main(args, config):
    """
    Main function of the program PyCGTOOL.

    Performs the complete AA->CG mapping and outputs a GROMACS forcefield directory.

    :param args: Arguments from argparse
    :param config: Configuration dictionary
    """
    frame = Frame(gro=args.gro, xtc=args.xtc)

    if args.bnd:
        bonds = BondSet(args.bnd)

    if args.map:
        mapping = Mapping(args.map)

    # Main loop - perform mapping and measurement on every frame in XTC
    while True:
        if args.map:
            cgframe = mapping.apply(frame, exclude={"SOL"})
        else:
            cgframe = frame

        if args.bnd:
            bonds.apply(cgframe)

        if not frame.next_frame():
            break

    if args.map:
        cgframe.output("out.gro", format=config["output"])

    if args.bnd:
        bonds.boltzmann_invert()
        for mol in bonds:
            print("Bonds in {0}:".format(mol))
            for bond in bonds[mol]:
                print(len(bond.values), bond.eqm, bond.fconst)

        bonds.write_itp("out.itp", mapping=mapping)
        ff = ForceField("fftest.ff")
        ff.write_rtp("test.rtp", mapping, bonds)


def map_only(args, config):
    """
    Perform AA->CG mapping and output coordinate file.

    :param args: Program arguments
    """
    frame = Frame(gro=args.gro)
    mapping = Mapping(args.map)
    cgframe = mapping.apply(frame, exclude={"SOL"})
    cgframe.output("out.gro", format=config["output"])


def interactive(default_config={}):
    """
    Read in options in interactive terminal mode.  Take defaults as options argument.

    :param default_config: Optional default configuration
    :return: Dictionary of configuration options
    """
    # TODO add interactive terminal mode to select options eg output format instead of CGTOOL config file
    config = {key: value for key, value in default_config.items()}
    while True:
        line = input(">>")
        if not line:
            break
        try:
            key, value = map(lambda x: x.lower(), line.split())
            if key in config:
                config[key] = value
            else:
                print("Invalid option")
        except ValueError:
            print("Must provide key value pair")
    return config


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform coarse-grain mapping of atomistic trajectory")
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-g', '--gro', type=str, required=True, help="GROMACS GRO file")
    required.add_argument('-x', '--xtc', type=str, help="GROMACS XTC file")
    required.add_argument('-m', '--map', type=str, help="Mapping file")
    required.add_argument('-b', '--bnd', type=str, help="Bonds file")

    parser.add_argument('-i', '--interactive', default=False, action='store_true')

    args = parser.parse_args()
    print("Using GRO: {0}".format(args.gro))
    print("Using XTC: {0}".format(args.xtc))

    default_config = {"output": "gro",
                      "map-only": "no",
                      "map-center": "geom"}
    if args.interactive:
        config = interactive(default_config)
    else:
        config = default_config

    if config["map-only"] == "yes":
        map_only(args, config)
    else:
        main(args, config)
