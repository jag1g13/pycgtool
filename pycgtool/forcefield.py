"""
This module contains a single class ForceField used to output a GROMACS .ff forcefield.
"""

import os
import shutil
import functools

from .util import dir_up
from .parsers import ITP


class ForceField:
    """
    Class used to output a GROMACS .ff forcefield
    """
    def __init__(self, name):
        """
        Open a named forcefield directory.  If it does not exist it is created.

        :param name: Forcefield name to open/create
        """
        self.dirname = "ff{0}.ff".format(name)
        os.makedirs(self.dirname, exist_ok=True)

        with open(os.path.join(self.dirname, "forcefield.itp"), "w") as itp:
            print("#define _FF_PYCGTOOL_{0}".format(name), file=itp)
            print('#include "martini_v2.2.itp"', file=itp)

        dist_dat_dir = os.path.join(dir_up(os.path.realpath(__file__), 2), "data")
        # Copy main MARTINI itp
        martini_itp = os.path.join(dist_dat_dir, "martini_v2.2.itp")
        shutil.copyfile(martini_itp, os.path.join(self.dirname, "martini_v2.2.itp"))
        # Copy water models
        shutil.copyfile(os.path.join(dist_dat_dir, "watermodels.dat"),
                        os.path.join(self.dirname, "watermodels.dat"))
        shutil.copyfile(os.path.join(dist_dat_dir, "w.itp"),
                        os.path.join(self.dirname, "w.itp"))

        # Create atomtypes.atp required for correct masses with pdb2gmx
        atomtypes_atp = os.path.join(self.dirname, "atomtypes.atp")
        with ITP(martini_itp) as itp, open(atomtypes_atp, 'w') as atomtypes:
            for toks in itp["atomtypes"]:
                print(" ".join(toks), file=atomtypes)

        with open(os.path.join(self.dirname, "forcefield.doc"), "w") as doc:
            print("PyCGTOOL produced MARTINI force field - {0}".format(name), file=doc)

    def write(self, filename, mapping, bonds):
        nterms, cterms, bothterms = self._write_rtp(filename, mapping, bonds)
        self._write_r2b(filename, nterms, cterms, bothterms)

    def _write_rtp(self, filename, mapping, bonds):
        """
        Write a GROMACS .rtp file.

        This file defines the residues present in the forcefield and allows pdb2gmx to be used.

        :param filename: Name of the .rtp file to create, N.B. .rtp is appended here
        :param mapping: AA->CG mapping from which to collect molecules
        :param bonds: BondSet from which to collect bonds
        """
        def write_bond_angle_dih(bonds, section_header, file, multiplicity=None):
            if bonds:
                print("  [ {0:s} ]".format(section_header), file=file)
            for bond in bonds:
                line = "    " + " ".join(["{0:>4s}".format(atom) for atom in bond.atoms])
                line += " {0:12.5f} {1:12.5f}".format(bond.eqm, bond.fconst)
                if multiplicity is not None:
                    line += " {0:4d}".format(multiplicity)
                print(line, file=file)

        def any_starts_with(iterable, char):
            """
            Return True if any atoms of any bonds in molecule start with 'char'.

            i.e. if char='-' or '+' is part of polymer.

            :param iterable: Iterable of bond entries to check
            :param char: Char to check each atom name for startswith, in '-+'
            :return: True if any atom name in molecule bonds starts with char, else False
            """
            recurse = functools.partial(any_starts_with, char=char)
            if type(iterable) is str:
                return iterable.startswith(char)
            else:
                return any(map(recurse, iterable))

        def write_residue(name, rtp, strip=None, prepend=""):
            print("[ {0} ]".format(prepend + name), file=rtp)

            print("  [ atoms ]", file=rtp)
            for bead in mapping[name]:
                #          name  type  charge  chg-group
                print("    {:>4s} {:>4s} {:3.6f} {:4d}".format(
                    bead.name, bead.type, bead.charge, 0
                ), file=rtp)

            needs_terminal_entry = [False, False]

            get_bond_functions = [("bonds", functools.partial(bonds.get_bond_lengths, with_constr=True)),
                                  ("angles", bonds.get_bond_angles),
                                  ("dihedrals", bonds.get_bond_dihedrals)]

            for get_bond in get_bond_functions:
                bond_tmp = get_bond[1](name)
                if strip is not None:
                    bond_tmp = [bond for bond in bond_tmp if not any_starts_with(bond, strip)]
                write_bond_angle_dih(bond_tmp, get_bond[0], rtp,
                                     multiplicity=1 if get_bond[0] == "dihedrals" else None)
                needs_terminal_entry[0] |= any_starts_with(bond_tmp, "-")
                needs_terminal_entry[1] |= any_starts_with(bond_tmp, "+")

            return needs_terminal_entry

        n_terms = set()
        c_terms = set()
        both_terms = set()

        with open(os.path.join(self.dirname, filename + ".rtp"), "w") as rtp:
            print("[ bondedtypes ]", file=rtp)
            print(("{:4d}" * 8).format(1, 1, 1, 1, 1, 1, 0, 0), file=rtp)

            for mol in mapping:
                # Skip molecules not listed in bonds
                if mol not in bonds:
                    continue

                needs_terminal_entry = write_residue(mol, rtp)
                if needs_terminal_entry[0]:
                    write_residue(mol, rtp, strip="-", prepend="N")
                    n_terms.add(mol)
                if needs_terminal_entry[1]:
                    write_residue(mol, rtp, strip="+", prepend="C")
                    c_terms.add(mol)
                if all(needs_terminal_entry):
                    write_residue(mol, rtp, strip=("-", "+"), prepend="2")
                    both_terms.add(mol)

        return n_terms, c_terms, both_terms

    def _write_r2b(self, filename, n_terms, c_terms, both_terms):
        with open(os.path.join(self.dirname, filename + ".r2b"), "w") as r2b:
            print("; rtp residue to rtp building block table", file=r2b)
            print(";     main  N-ter C-ter 2-ter", file=r2b)

            for resname in set.union(n_terms, c_terms, both_terms):
                nter_str = ("N" + resname) if resname in n_terms else "-"
                cter_str = ("C" + resname) if resname in c_terms else "-"
                both_ter_str = ("2" + resname) if resname in both_terms else "-"
                print("{0:5s} {0:5s} {1:5s} {2:5s} {3:5s}".format(resname, nter_str, cter_str, both_ter_str), file=r2b)
