"""
This module contains classes required to perform an atomistic to coarse-grain mapping.

The Mapping class contains a dictionary of lists of BeadMaps.
Each list corresponds to a single molecule.
"""

import numpy as np

from .frame import Atom, Residue, Frame
from .parsers.cfg import CFG
from .util import dist_with_pbc

try:
    from numba import jit
except ImportError:
    from .util import jit_dummy as jit

np.seterr(all="raise")


class BeadMap(Atom):
    """
    POD class holding values relating to the AA->CG transformation for a single bead.
    """
    __slots__ = ["name", "type", "atoms", "charge", "mass", "weights"]

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
        Atom.__init__(self, name=name, type=type, charge=charge, mass=mass)
        self.atoms = atoms
        self.weights = {"geom": None,
                        "mass": np.ones(len(self.atoms), dtype=np.float32)}

    def __iter__(self):
        """
        Iterate through the atom names from which the bead is made up.

        :return: Iterator over atoms
        """
        return iter(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self, item):
        return self.atoms[item]


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
    def __init__(self, filename, options, itp=None):
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
                    assert atoms, "Bead {0} specification contains no atoms".format(name)
                    newbead = BeadMap(name=name, type=typ, atoms=atoms)
                    molmap.append(newbead)

        # TODO this only works with one moleculetype in one itp - extend this
        if itp is not None:
            with CFG(itp) as itp:
                atoms = {}
                for toks in itp["atoms"]:
                    # Store charge and mass
                    atoms[toks[4]] = (float(toks[6]), float(toks[7]))

                for bead in self._mappings[itp["moleculetype"][0][0]]:
                    bead.weights["mass"] = np.array([[atoms[atom][1]] for atom in bead], dtype=np.float32)
                    for atom in bead:
                        bead.charge += atoms[atom][0]
                        bead.mass += atoms[atom][1]

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
            if aares.name not in self._mappings:
                continue
            if exclude is not None and aares.name in exclude:
                continue

            res = self._apply_res_pbc(aares, frame.box)

            cgframe.natoms += len(res)
            cgframe.residues.append(res)
        return cgframe

    def _apply_res_pbc(self, aares, box):
        """
        Apply mapping transformation to a single residue to allow multithreading.

        :param aares: Atomistic residue to apply mapping
        :param box: Cubic periodic box vectors
        :return: A single coarse grained residue
        """

        molmap = self._mappings[aares.name]
        res = Residue(name=aares.name, num=aares.num)
        res.atoms = [Atom(name=bead.name, type=bead.type, charge=bead.charge, mass=bead.mass) for bead in molmap]

        # Perform mapping
        for i, (bead, bmap) in enumerate(zip(res, molmap)):
            res.name_to_num[bead.name] = i
            bead.charge = bmap.charge
            bead.mass = bmap.mass

            ref_coords = aares[bmap[0]].coords
            coords = np.array([aares[atom].coords for atom in bmap], dtype=np.float32)

            weights = bmap.weights[self._map_center]
            if weights is None:
                bead.coords = calc_coords(ref_coords, coords, box)
            else:
                bead.coords = calc_coords_weight(ref_coords, coords, box, weights)

        return res


@jit
def calc_coords_weight(ref_coords, coords, box, weights):
    coords = dist_with_pbc(ref_coords, coords, box)
    coords = np.sum(weights * coords, axis=0)
    coords /= np.sum(weights)
    coords += ref_coords
    return coords


@jit
def calc_coords(ref_coords, coords, box):
    coords = dist_with_pbc(ref_coords, coords, box)
    n = len(coords)
    coords = np.sum(coords, axis=0)
    coords /= n
    coords += ref_coords
    return coords

