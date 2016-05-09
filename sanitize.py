#!/usr/bin/env python3
import argparse


def main(args):
    with open(args.i) as infile, open(args.o, "w") as outfile:
        section = ""
        expected_toks = {"bonds":       5,
                         "angles":      6,
                         "dihedrals":   8,
                         "constraints": 4,
                         "pairs":      -1}
        has_constr = False
        constr_lines = []

        for line in infile:
            line = line.strip("\n")
            if line.startswith("["):
                section = line.strip("[ ]")
                if args.fc is not None and section == "constraints":
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

            elif line.startswith("#") or line.startswith(";") or line == "":
                print(line, file=outfile)
                continue

            elif len(line.split()) < expected_toks.get(section, 0):
                continue

            elif expected_toks.get(section, 0) < 0:
                continue

            # Convert high force constants to constraints
            if args.fc is not None and section == "bonds" and float(line.split()[4]) >= args.fc:
                constr_lines.append(line[:line.rindex(line.split()[4]) - 1])
                continue

            print(line, file=outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sanitize PDB2GMX output topology.")
    parser.add_argument('-i', type=str, required=True, help="Input topology file")
    parser.add_argument('-o', type=str, required=True, help="Output topology file")
    parser.add_argument('-fc', type=float, default=None, help="Convert high force constants to constraints")

    args = parser.parse_args()

    main(args)
