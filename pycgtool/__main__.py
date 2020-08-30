#!/usr/bin/env python3

import argparse
import sys

try:
    from pycgtool.pycgtool import main, map_only
    from pycgtool.interface import Options
    from pycgtool.functionalforms import FunctionalForms
except SyntaxError:
    raise RuntimeError("PyCGTOOL requires Python 3.2 or greater")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform coarse-grain mapping of atomistic trajectory")
    input_files = parser.add_argument_group("Input files")
    input_files.add_argument('-g', '--gro', type=str, required=True, help="GROMACS GRO file")
    input_files.add_argument('-m', '--map', type=str, help="Mapping file")
    input_files.add_argument('-x', '--xtc', type=str, help="GROMACS XTC file")
    input_files.add_argument('-b', '--bnd', type=str, help="Bonds file")
    input_files.add_argument('-i', '--itp', type=str, help="GROMACS ITP file")

    parser.add_argument('--advanced', default=False, action='store_true', help="Show advanced options menu")
    parser.add_argument('--outputxtc', default=False, action='store_true', help="Output a pseudo-CG trajectory")
    parser.add_argument('--quiet', default=False, action='store_true', help="Hide progress bars")
    input_files.add_argument('--begin', type=int, default=0, help="Frame number to begin")
    input_files.add_argument('--end', type=int, default=-1, help="Frame number to end")

    func_forms = FunctionalForms()

    args = parser.parse_args()
    config = Options([
        ("output_name", "out"),
        ("output", "gro"),
        ("output_xtc", args.outputxtc),
        ("map_only", not bool(args.bnd)),
        ("map_center", "geom"),
        ("constr_threshold", 100000),
        ("dump_measurements", bool(args.bnd) and not bool(args.map)),
        ("dump_n_values", 10000),
        ("output_forcefield", False),
        ("temperature", 310),
        ("default_fc", False),
        ("generate_angles", True),
        ("generate_dihedrals", False),
        ("length_form", "harmonic"),
        ("angle_form", "cosharmonic"),
        ("dihedral_form", "harmonic")
    ], args)

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
        main(args, config)
