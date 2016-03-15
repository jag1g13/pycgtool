"""
This module contains classes required to perform an atomistic to coarse-grain mapping.

The Mapping class contains a dictionary of lists of BeadMaps.
Each list corresponds to a single molecule.
"""

import numpy as np

from .frame import Atom, Residue, Frame
from .parsers.cfg import CFG

np.seterr(all="raise")


class BeadMap:
    """
    POD class holding values relating to the AA->CG transformation for a single bead.
    """
    __slots__ = ["name", "typ", "type", "atoms", "charge", "mass"]

    def __init__(self, name=None, type=None, atoms=None, charge=0, mass=0):
        """
        Create a single bead mapping

        :param name: The name of the bead
        :param type: The bead type
        :param atoms: The atom names from which the bead is made up
        :param charge: The net charge on the bead
        :param mass: The total bead mass
        :return: Instance of BeadMap
        """
        self.name = name
        self.type = type
        self.atoms = atoms
        self.charge = charge
        self.mass = mass

    def __iter__(self):
        """
        Iterate through the atom names from which the bead is made up.

        :return: Iterator over atoms
        """
        return iter(self.atoms)

    def __len__(self):
        return len(self.atoms)


class EmptyBeadError(Exception):
    """
    Exception used to indicate that none of the required atoms are present.
    """
    pass


def coordinate_weight(center, atom):
    centers = {"geom": lambda at: at.coords,
               "mass": lambda at: at.coords * at.mass}
    try:
        return centers[center](atom)
    except KeyError:
        raise ValueError("Invalid map-center type '{0}'".format(center))


class Mapping:
    """
    Class used to perform the AA->CG mapping.

    Contains a dictionary of lists of BeadMaps.  Each list corresponds to a single molecule.
    """
    def __init__(self, filename, options):
        """
        Read in the AA->CG mapping from a file.

        :param filename: File from which to read mapping
        :return: Instance of Mapping
        """
        self._mappings = {}

        self._map_center = options.map_center

        with CFG(filename) as cfg:
            for mol in cfg:
                self._mappings[mol.name] = []
                molmap = self._mappings[mol.name]
                for name, typ, *atoms in mol:
                    newbead = BeadMap(name=name, type=typ, atoms=atoms)
                    molmap.append(newbead)

    def __len__(self):
        return len(self._mappings)

    def __contains__(self, item):
        return item in self._mappings

    def __getitem__(self, item):
        return self._mappings[item]

    def __iter__(self):
        return iter(self._mappings)

    def apply(self, frame, exclude=None):
        """
        Apply the AA->CG mapping to an atomistic Frame.

        :param frame: Frame to which mapping will be applied
        :param exclude: Set of molecule names to exclude from mapping - e.g. solvent
        :return: A new Frame instance containing the CG frame
        """
        cgframe = Frame()
        cgframe.name = frame.name
        cgframe.box = frame.box
        cgframe.natoms = 0
        cgframe.residues = []

        for aares in frame:
            if exclude is not None and aares.name in exclude:
                continue
            try:
                molmap = self._mappings[aares.name]
            except KeyError:
                # Residue not in mapping, assume user doesn't want to map it (e.g. solvent)
                continue
            res = Residue(name=aares.name, num=aares.num)
            res.atoms = [Atom(name=atom.name, type=atom.type) for atom in molmap]

            # Perform mapping
            for i, (bead, bmap) in enumerate(zip(res, molmap)):
                res.name_to_num[bead.name] = i
                bead.coords = np.zeros(3)
                n = 0
                for atom in bmap:
                    try:
                        bead.coords += coordinate_weight(self._map_center, aares[atom])
                        n += 1
                    except KeyError:
                        # Atom not in residue, happens in terminal residues of polymer
                        pass

                try:
                    bead.coords /= n
                except FloatingPointError:
                    # Bead contains no atoms, happens in terminal residue of polymer
                    pass
                cgframe.natoms += 1

            cgframe.residues.append(res)
        return cgframe
