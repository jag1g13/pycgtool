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
        for line in infile:
            line = line.strip("\n")
            if line.startswith("["):
                section = line.strip("[ ]")
            elif line.startswith("#") or line.startswith(";"):
                pass
            elif len(line.split()) < expected_toks.get(section, 0):
                continue

            if expected_toks.get(section, 0) < 0:
                continue

            print(line, file=outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Sanitize PDB2GMX output topology.")
    parser.add_argument('-i', type=str, required=True, help="Input topology file")
    parser.add_argument('-o', type=str, required=True, help="Output topology file")

    args = parser.parse_args()

    main(args)
