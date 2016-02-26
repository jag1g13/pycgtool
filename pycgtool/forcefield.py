"""
This module contains a single class ForceField used to output a GROMACS .ff forcefield.
"""

import os

# Python 3.2 doesn't have FileExistsError
try:
    raise FileExistsError
except NameError:
    class FileExistsError(Exception):
        pass
except FileExistsError:
    pass


class ForceField:
    """
    Class used to output a GROMACS .ff forcefield
    """
    def __init__(self, name):
        """
        Open a forcefield directory at name.  If it does not exist it is created.

        :param name: Directory name to open/create
        """
        try:
            os.makedirs(name)
        except FileExistsError as e:
            if not os.path.isdir(name):
                raise e
        except OSError as e:
            if e.errno != 17 or not os.path.isdir(name):
                raise e

        self.dirname = name
        with open(os.path.join(self.dirname, "forcefield.itp"), "w") as itp:
            print("#define _FF_PYCGTOOL", file=itp)
            print("\n[ defaults ]", file=itp)
            print("{:4d} {:4d} {:4s} {:2.4f} {:2.4f}".format(
                1, 1, "no", 1.0, 1.0
            ), file=itp)

        with open(os.path.join(self.dirname, "forcefield.doc"), "w") as doc:
            print("PyCGTOOL produced force field", file=doc)

    def write_atp(self, mapping):
        """
        Write a GROMACS forcefield .atp file.

        This file lists the atomtypes used in the forcefield.

        :param mapping: AA->CG mapping from which to collect atomtypes
        """
        with open(os.path.join(self.dirname, "atomtypes.atp"), "w") as atp:
            types = set()
            for mol in mapping:
                for bead in mapping[mol]:
                    if bead.type not in types:
                        print("{:4s} {:8.3f}".format(bead.type, bead.mass), file=atp)
                        types.add(bead.type)

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
            print(("{:4d}"*4).format(1, 1, 1, 1), file=rtp)

            for mol in mapping:
                print("[ {0} ]".format(mol), file=rtp)

                print(" [ atoms ]", file=rtp)
                for i, bead in enumerate(mapping[mol]):
                    #        name  type  charge  c-group
                    print("  {:4s} {:4s} {:3.6f} {:4d}".format(
                        bead.name, bead.type, bead.charge, i+1
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
