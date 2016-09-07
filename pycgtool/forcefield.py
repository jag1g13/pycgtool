"""
This module contains a single class ForceField used to output a GROMACS .ff forcefield.
"""

import os
import shutil
import functools

from .util import dir_up
from .parsers.cfg import CFG


class ForceField:
    """
    Class used to output a GROMACS .ff forcefield
    """
    # TODO output r2b file for ends of chain
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
        shutil.copyfile(os.path.join(dist_dat_dir, "martini_v2.2.itp"),
                        os.path.join(self.dirname, "martini_v2.2.itp"))
        # Copy water models
        shutil.copyfile(os.path.join(dist_dat_dir, "watermodels.dat"),
                        os.path.join(self.dirname, "watermodels.dat"))
        shutil.copyfile(os.path.join(dist_dat_dir, "w.itp"),
                        os.path.join(self.dirname, "w.itp"))

        # Create atomtypes.atp required for correct masses with pdb2gmx
        with CFG(os.path.join(dist_dat_dir, "martini_v2.2.itp"), allow_duplicate=True) as cfg,\
                open(os.path.join(self.dirname, "atomtypes.atp"), 'w') as atomtypes:
            for toks in cfg["atomtypes"]:
                print(" ".join(toks), file=atomtypes)

        with open(os.path.join(self.dirname, "forcefield.doc"), "w") as doc:
            print("PyCGTOOL produced MARTINI force field - {0}".format(name), file=doc)

    def write_rtp(self, name, mapping, bonds):
        """
        Write a GROMACS .rtp file.

        This file defines the residues present in the forcefield and allows pdb2gmx to be used.

        :param name: Name of the .rtp file to create
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

        def strip_polymer_bonds(bonds, char):
            return [bond for bond in bonds if not any_starts_with(bond, char)]

        # def write_residue(name):

        with open(os.path.join(self.dirname, name), "w") as rtp:
            print("[ bondedtypes ]", file=rtp)
            print(("{:4d}" * 8).format(1, 1, 1, 1, 1, 1, 0, 0), file=rtp)

            for mol in mapping:
                # Skip molecules without bonds
                if mol not in bonds:
                    continue

                print("[ {0} ]".format(mol), file=rtp)

                print("  [ atoms ]", file=rtp)
                for bead in mapping[mol]:
                    #          name  type  charge  chg-group
                    print("    {:>4s} {:>4s} {:3.6f} {:4d}".format(
                        bead.name, bead.type, bead.charge, 0
                    ), file=rtp)

                needs_terminal_entry = [False, False]

                bond_tmp = bonds.get_bond_lengths(mol, with_constr=True)
                write_bond_angle_dih(bond_tmp, "bonds", rtp)
                needs_terminal_entry[0] |= any_starts_with(bond_tmp, "-")
                needs_terminal_entry[1] |= any_starts_with(bond_tmp, "+")

                bond_tmp = bonds.get_bond_angles(mol)
                write_bond_angle_dih(bond_tmp, "angles", rtp)
                needs_terminal_entry[0] |= any_starts_with(bond_tmp, "-")
                needs_terminal_entry[1] |= any_starts_with(bond_tmp, "+")

                bond_tmp = bonds.get_bond_dihedrals(mol)
                write_bond_angle_dih(bond_tmp, "dihedrals", rtp, multiplicity=1)
                needs_terminal_entry[0] |= any_starts_with(bond_tmp, "-")
                needs_terminal_entry[1] |= any_starts_with(bond_tmp, "+")

                # print(mol, needs_terminal_entry)
                # bond_tmp = bonds.get_bond_lengths(mol, with_constr=True)
                # if needs_terminal_entry[0]:
                #     print(strip_polymer_bonds(bond_tmp, "-"))

