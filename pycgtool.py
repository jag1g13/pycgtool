#!/usr/bin/env python3

import argparse

from pycgtool.frame import Frame
from pycgtool.mapping import Mapping
from pycgtool.bondset import BondSet
from pycgtool.forcefield import ForceField


def main(args):
    frame = Frame(gro=args.gro, xtc=args.xtc)

    if args.bnd is not None:
        bonds = BondSet(args.bnd)

    if args.map is not None:
        mapping = Mapping(args.map)

    # Main loop - perform mapping and measurement on every frame in XTC
    while True:
        if args.map is None:
            cgframe = frame
        else:
            cgframe = mapping.apply(frame, exclude={"SOL"})

        if args.bnd is not None:
            bonds.apply(cgframe)

        if not frame.next_frame():
            break

    if args.map is not None:
        cgframe.output_gro("out.gro")

    if args.bnd is not None:
        bonds.boltzmann_invert()
        for mol in bonds:
            print("Bonds in {0}:".format(mol))
            for bond in bonds[mol]:
                print(len(bond.values), bond.eqm, bond.fconst)

        bonds.write_itp("out.itp", mapping=mapping)
        ff = ForceField("fftest.ff")
        ff.write_rtp("test.rtp", mapping, bonds)
        ff.write_atp(mapping)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform coarse-grain mapping of atomistic trajectory")
    required = parser.add_argument_group("Required arguments")
    required.add_argument('-g', '--gro', type=str, required=True, help="GROMACS GRO file")
    required.add_argument('-x', '--xtc', type=str, help="GROMACS XTC file")
    required.add_argument('-m', '--map', type=str, help="Mapping file")
    required.add_argument('-b', '--bnd', type=str, help="Bonds file")

    args = parser.parse_args()
    print("Using GRO: {0}".format(args.gro))
    print("Using XTC: {0}".format(args.xtc))

    main(args)
