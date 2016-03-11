"""
This module contains a single class ForceField used to output a GROMACS .ff forcefield.
"""

import os
import shutil

from .util import dir_up
from .parsers.cfg import CFG


class ForceField:
    """
    Class used to output a GROMACS .ff forcefield
    """
    def __init__(self, name):
        """
        Open a forcefield directory at name.  If it does not exist it is created.

        :param name: Directory name to open/create
        """
        os.makedirs(name, exist_ok=True)

        self.dirname = name
        with open(os.path.join(self.dirname, "forcefield.itp"), "w") as itp:
            print("#define _FF_PYCGTOOL", file=itp)
            print('#include "martini_v2.2.itp"', file=itp)

        # Copy main MARTINI itp
        shutil.copyfile(os.path.join(dir_up(os.path.realpath(__file__), 2), "data", "martini_v2.2.itp"),
                        os.path.join(self.dirname, "martini_v2.2.itp"))
        # Copy water models
        shutil.copyfile(os.path.join(dir_up(os.path.realpath(__file__), 2), "data", "watermodels.dat"),
                        os.path.join(self.dirname, "watermodels.dat"))
        shutil.copyfile(os.path.join(dir_up(os.path.realpath(__file__), 2), "data", "w.itp"),
                        os.path.join(self.dirname, "w.itp"))

        # Create atomtypes.atp required for correct masses with pdb2gmx
        with CFG(os.path.join(dir_up(os.path.realpath(__file__), 2), "data", "martini_v2.2.itp"), allow_duplicate=True) as cfg,\
                open(os.path.join(self.dirname, "atomtypes.atp"), 'w') as atomtypes:
            for toks in cfg["atomtypes"]:
                print(" ".join(toks), file=atomtypes)

        with open(os.path.join(self.dirname, "forcefield.doc"), "w") as doc:
            print("PyCGTOOL produced MARTINI force field", file=doc)

    def write_rtp(self, name, mapping, bonds):
        """
        Write a GROMACS .rtp file.

        This file defines the residues present in the forcefield and allows pdb2gmx to be used.

        :param name: Name of the .rtp file to create
        :param mapping: AA->CG mapping from which to collect molecules
        :param bonds: BondSet from which to collect bonds
        """
        # TODO print everything to file
        with open(os.path.join(self.dirname, name), "w") as rtp:
            print("[ bondedtypes ]", file=rtp)
            print(("{:4d}"*8).format(1, 1, 1, 1, 1, 1, 0, 0), file=rtp)

            for mol in mapping:
                try:
                    bonds[mol]
                except KeyError:
                    # Skip molecules without bonds
                    continue

                print("[ {0} ]".format(mol), file=rtp)

                print(" [ atoms ]", file=rtp)
                for bead in mapping[mol]:
                    #        name  type  charge  c-group
                    print("  {:4s} {:4s} {:3.6f} {:4d}".format(
                        bead.name, bead.type, bead.charge, 0
                    ), file=rtp)

                bond_tmp = [bond for bond in bonds[mol] if len(bond) == 2]
                if bond_tmp:
                    print(" [ bonds ]", file=rtp)
                for bond in bond_tmp:
                    print("  {:4s} {:4s} {:8.3f} {:8.3f}".format(
                        bond.atoms[0], bond.atoms[1],
                        bond.eqm, bond.fconst
                    ), file=rtp)

                bond_tmp = [bond for bond in bonds[mol] if len(bond) == 3]
                if bond_tmp:
                    print(" [ angles ]", file=rtp)
                for bond in bond_tmp:
                    print("  {:4s} {:4s} {:4s} {:8.3f} {:8.3f}".format(
                        bond.atoms[0], bond.atoms[1], bond.atoms[2],
                        bond.eqm, bond.fconst
                    ), file=rtp)

                bond_tmp = [bond for bond in bonds[mol] if len(bond) == 4]
                if bond_tmp:
                    print(" [ dihedrals ]", file=rtp)
                for bond in bond_tmp:
                    print("  {:4s} {:4s} {:4s} {:4s} {:8.3f} {:8.3f} {:4d}".format(
                        bond.atoms[0], bond.atoms[1], bond.atoms[2], bond.atoms[3],
                        bond.eqm, bond.fconst, 1
                    ), file=rtp)
