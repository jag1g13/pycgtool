#!/usr/bin/env python3

import argparse

from pycgtool.frame import Frame, Mapping, Measure


def main(args):
    frame = Frame(gro=args.gro, xtc=args.xtc)

    if args.bnd is not None:
        bonds = Measure(args.bnd)

    if args.map is not None:
        mapping = Mapping(args.map)
        cgframe = mapping.apply(frame, exclude={"SOL"})

    # Main loop - perform mapping and measurement on every frame in XTC
    while frame.next_frame():
        if args.map is None:
            cgframe = frame
        else:
            cgframe = mapping.apply(frame, exclude={"SOL"})

        if args.bnd is not None:
            bonds.apply(cgframe)

    if args.bnd is not None:
        for mol in bonds:
            print("Bonds in {0}:".format(mol))
            for bond in bonds[mol]:
                # print(len(bond.values))
                try:
                    print(sum(bond.values) / len(bond.values))
                except ZeroDivisionError:
                    print("Bond has no values")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform coarse-grain mapping of atomistic trajectory")
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-g', '--gro', type=str, required=True, help="GROMACS GRO file")
    required.add_argument('-x', '--xtc', type=str, required=True, help="GROMACS XTC file")
    required.add_argument('-m', '--map', type=str, help="Mapping file")
    required.add_argument('-b', '--bnd', type=str, help="Bonds file")

    args = parser.parse_args()
    print("Using GRO: {0}".format(args.gro))
    print("Using XTC: {0}".format(args.xtc))

    main(args)
