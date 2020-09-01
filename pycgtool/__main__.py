#!/usr/bin/env python3

import argparse
import logging
import pathlib

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
                             help="Frame number to begin")
    input_files.add_argument('--end', type=int, default=-1,
                             help="Frame number to end")

    parser.add_argument('--out-dir', default='.', type=str,
                        help="Directory where output files should be placed")
    parser.add_argument('--outputxtc', default=False, action='store_true',
                        help="Output a pseudo-CG trajectory")
    parser.add_argument('--quiet', default=False, action='store_true',
                        help="Hide progress bars")

    # Advanced options
    advanced_arguments = parser.add_argument_group("Advanced options")

    advanced_arguments.add_argument("--output_name", default="out",
                                    help="Base name of output files")
    advanced_arguments.add_argument("--output", default="gro",
                                    help="Coordinate output format")
    advanced_arguments.add_argument("--map_only", default=False, action="store_true",
                                    help="Run in mapping-only mode")
    advanced_arguments.add_argument("--map_center", default="geom",
                                    choices=["geom", "mass", "first"],
                                    help="Mapping method")
    advanced_arguments.add_argument("--virtual_map_center", default="geom",
                                    choices=["geom", "mass"],
                                    help="Virtual site mapping method")
    advanced_arguments.add_argument("--constr_threshold", type=float, default=100000,
                                    help="Convert stiff bonds to contraints over [value]")
    advanced_arguments.add_argument("--dump_measurements", default=False, action="store_true",
                                    help="Whether to output bond measurements")
    advanced_arguments.add_argument("--dump_n_values", type=int, default=10000,
                                    help="How many measurements to output")
    advanced_arguments.add_argument("--output_forcefield", default=False, action='store_true',
                                    help="Output GROMACS forefield directory?")
    advanced_arguments.add_argument("--temperature", type=float, default=310,
                                    help="Temperature of reference simulation")
    advanced_arguments.add_argument("--default_fc", default=False, action='store_true',
                                    help="Use default MARTINI force constants?")
    advanced_arguments.add_argument("--generate_angles", default=True, action='store_false',
                                    help="Generate angles from bonds")
    advanced_arguments.add_argument("--generate_dihedrals", default=False, action="store_true",
                                    help="Generate dihedrals from bonds")
    advanced_arguments.add_argument("--length_form", default="harmonic",
                                    help="Form of bond potential")
    advanced_arguments.add_argument("--angle_form", default="cosharmonic",
                                    help="Form of angle potential")
    advanced_arguments.add_argument("--dihedral_form", default="harmonic",
                                    help="Form of dihedral potential")
    # yapf: enable

    # Populate functional forms dictionary - ignore returned value
    _ = FunctionalForms()

    args = parser.parse_args()

    if not args.dump_measurements:
        args.dump_measurements = bool(args.bnd) and not bool(args.map)
    if not args.map_only:
        args.map_only = not bool(args.bnd)

    config = Options([
        ("output_name", args.output_name),
        ("output", args.output),
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
    ], args)  # yapf: disable

    if not args.map and not args.bnd:
        parser.error("One or both of -m and -b is required.")

    print("Using GRO: {0}".format(args.gro))
    print("Using XTC: {0}".format(args.xtc))

    if config.map_only:
        map_only(args, config)

    else:
        full_run(args, config)


if __name__ == "__main__":
    main()
