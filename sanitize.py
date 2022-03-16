#!/usr/bin/env python3
import argparse

EXPECTED_TOKS = {
    "bonds": 5,
    "angles": 6,
    "dihedrals": 8,
    "constraints": 4,
    "pairs": -1
}


def main(args):
    with open(args.input) as infile, open(args.output, "w") as outfile:
        section = ""
        has_constr = False
        constr_lines = []

        for line in infile:
            line = line.strip("\n")
            if line.startswith("["):
                section = line.strip("[ ]")
                if args.constraint_threshold is not None and section == "constraints":
                    print(line, file=outfile)
                    has_constr = True
                    for c_line in constr_lines:
                        print(c_line, file=outfile)
                    continue

                if section == "system" and not has_constr and constr_lines:
                    print("[ constraints ]", file=outfile)
                    for c_line in constr_lines:
                        print(c_line, file=outfile)
                    print(file=outfile)
                    has_constr = True

                print(line, file=outfile)
                continue

            if line.startswith("#") or line.startswith(";") or line == "":
                print(line, file=outfile)
                continue

            if len(line.split()) < EXPECTED_TOKS.get(section, 0):
                continue

            if EXPECTED_TOKS.get(section, 0) < 0:
                continue

            # Convert high force constants to constraints
            if (args.constraint_threshold is not None and section == "bonds"
                    and float(line.split()[4]) >= args.constraint_threshold):
                constr_lines.append(line[:line.rindex(line.split()[4]) - 1])
                continue

            print(line, file=outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Sanitize PDB2GMX output topology.")

    parser.add_argument('input', type=str, help="Input topology file")
    parser.add_argument('output', type=str, help="Output topology file")
    parser.add_argument('-c',
                        '--constraint-threshold',
                        type=float,
                        default=None,
                        help="Convert high force constants to constraints")

    main(parser.parse_args())
