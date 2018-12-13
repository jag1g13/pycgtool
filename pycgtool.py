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
    parser = argparse.ArgumentParser(description="Perform coarse-grain mapping of atomistic trajectory",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    advanced_arguments = parser.add_argument_group("Advanced options")
    input_files = parser.add_argument_group("Input files")

    input_files.add_argument('-g', '--gro', type=str, required=True, help="GROMACS GRO file")
    input_files.add_argument('-m', '--map', type=str, help="Mapping file")
    input_files.add_argument('-x', '--xtc', type=str, help="GROMACS XTC file")
    input_files.add_argument('-b', '--bnd', type=str, help="Bonds file")
    input_files.add_argument('-i', '--itp', type=str, help="GROMACS ITP file")

    parser.add_argument('--outputxtc', default=False, action='store_true', help="Output a pseudo-CG trajectory")
    parser.add_argument('--quiet', default=False, action='store_true', help="Hide progress bars")
    input_files.add_argument('--begin', type=int, default=0, help="Frame number to begin")
    input_files.add_argument('--end', type=int, default=-1, help="Frame number to end")
    #advanced options
    advanced_arguments.add_argument("--map_center", default="geom", choices=["geom", "mass", "first"], help="Mapping method")
    advanced_arguments.add_argument("--virtual_map_center", default="geom", choices=["geom", "mass"],
                             help="Virtual site mapping method")
    advanced_arguments.add_argument("--constr_threshold", type=float, default=100000,
                             help="Convert stiff bonds to contraints over [value]")
    advanced_arguments.add_argument("--dump_measurements", default=False, action="store_true",
                             help="Whether to output bond measurements")
    advanced_arguments.add_argument("--dump_n_values", type=int, default=10000, help="How many measurements to output")
    advanced_arguments.add_argument("--output_forcefield", default=False, action='store_true',
                             help="Output GROMACS forefield directory?")
    advanced_arguments.add_argument("--temperature", type=float, default=310, help="Temperature of reference simulation")
    advanced_arguments.add_argument("--default_fc", default=False, action='store_true', help="Use default MARTINI force constants?")
    advanced_arguments.add_argument("--generate_angles", default=True, action='store_false', help="Generate angles from bonds")
    advanced_arguments.add_argument("--generate_dihedrals", default=False, action="store_true", help="Generate dihedrals from bonds")
    advanced_arguments.add_argument("--length_form", default="harmonic", help="Form of bond potential")
    advanced_arguments.add_argument("--angle_form", default="cosharmonic", help="Form of angle potential")
    advanced_arguments.add_argument("--dihedral_form", default="harmonic", help="Form of dihedral potential")

    func_forms = FunctionalForms()

    args = parser.parse_args()
    if not args.dump_measurements:
        args.dump_measurements = bool(args.bnd) and not bool(args.map)
    if not args.map_only:
        args.map_only = not bool(args.bnd)

    config = Options([
        ("output_name", "out"),
        ("output", "gro"),
        ("output_xtc", args.outputxtc),
        ("map_only", args.map_only),
        ("map_center", args.map_center),
        ("virtual_map_center", args.virtual_map_center),
        ("constr_threshold", args.constr_threshold),
        ("dump_measurements", args.dump_measurements),
        ("dump_n_values", args.dump_n_values),
        ("output_forcefield", args.output_forcefield),
        ("temperature", args.temperature),
        ("default_fc", args.default_fc),
        ("generate_angles", args.generate_angles),
        ("generate_dihedrals", args.generate_dihedrals),
        ("length_form", args.length_form),
        ("angle_form", args.angle_form),
        ("dihedral_form", args.dihedral_form)
    ], args)

    if not args.map and not args.bnd:
        parser.error("One or both of -m and -b is required.")


    print("Using GRO: {0}".format(args.gro))
    print("Using XTC: {0}".format(args.xtc))

    if config.map_only:
        map_only(args, config)
    else:
        main(args, config)
