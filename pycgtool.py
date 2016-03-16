#!/usr/bin/env python3

import argparse
import curses
import curses.textpad

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


def interactive(config, stdscr):
    """
    Read in options in interactive terminal mode using curses.  Take defaults as options argument.

    :param config: Default configuration, will be mutated
    :param stdscr: Curses window to use as interface
    """
    stdscr.clear()
    stdscr.addstr(1, 1, "Using GRO: {0}".format(args.gro))
    stdscr.addstr(2, 1, "Using XTC: {0}".format(args.xtc))
    stdscr.box()
    stdscr.refresh()

    nrows = len(config)

    errscr = stdscr.derwin(2, curses.COLS - 3, nrows + 8, 1)
    errscr.border()

    window_config = stdscr.derwin(nrows + 2, curses.COLS - 3, 5, 1)
    window_config.box()
    window_config.refresh()
    window_keys = window_config.derwin(nrows, 20, 1, 0)
    window_config.vline(1, 18, curses.ACS_VLINE, nrows)
    window_vals = window_config.derwin(nrows, curses.COLS - 24, 1, 20)
    text_edit_wins = []
    text_inputs = []

    for i, (key, value) in enumerate(config):
        window_keys.addstr(i, 0, key)
        text_edit_wins.append(window_vals.derwin(1, 30, i, 0))
        text_edit_wins[-1].addstr(0, 0, str(value))
        text_inputs.append(curses.textpad.Textbox(text_edit_wins[-1]))

    stdscr.refresh()
    window_keys.refresh()
    for window in text_edit_wins:
        window.refresh()

    pos = 0
    nrows = len(config)
    move = {"KEY_UP": lambda x: max(x - 1, 0),
            "KEY_DOWN": lambda x: min(x + 1, nrows - 1),
            "KEY_LEFT": lambda x: x,
            "KEY_RIGHT": lambda x: x}

    while True:
        key = text_edit_wins[pos].getkey(0, 0)
        errscr.erase()
        if key in move:
            pos = move[key](pos)
        if key == "\n":
            if type(config[pos]) is bool:
                config.toggle_boolean_by_num(pos)
            else:
                val = text_inputs[pos].edit().strip()
                try:
                    config.set_by_num(pos, val)
                except ValueError:
                    errscr.addstr(0, 0, "Invalid value '{0}' for option".format(val))
                    errscr.addstr(1, 0, "Value has been reset".format(val))

            text_edit_wins[pos].erase()
            text_edit_wins[pos].addstr(0, 0, str(config[pos]))
            text_edit_wins[pos].refresh()

        errscr.refresh()
        if key == "q":
            break


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
    config = Options([("output", "gro"),
                      ("map_only", False),
                      ("map_center", "geom"),
                      ("constr_threshold", 100000),
                      ("dump_measurements", False),
                      ("output_forcefield", False)])
    if args.bnd is None:
        config.set("map_only", True)

    if args.interactive:
        def inter(stdscr):
            interactive(config, stdscr)
        curses.wrapper(inter)
    else:
        print("Using GRO: {0}".format(args.gro))
        print("Using XTC: {0}".format(args.xtc))

    if config.map_only:
        map_only(args, config)
    else:
        main(args, config)
