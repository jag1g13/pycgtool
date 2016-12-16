"""
This module contains classes required to perform an atomistic to coarse-grain mapping.

The Mapping class contains a dictionary of lists of BeadMaps.
Each list corresponds to a single molecule.
"""

import numpy as np
import logging
import json
import os

from .frame import Atom, Residue, Frame
from .parsers.cfg import CFG
from .util import dir_up

try:
    import numba
except ImportError:
    from .util import NumbaDummy
    numba = NumbaDummy()

np.seterr(all="raise")

logger = logging.getLogger(__name__)


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
        # NB: Mass weights are added in Mapping.__init__ if an itp file is provided
        self.weights = {"geom": np.array([[1. / len(atoms)] for _ in atoms], dtype=np.float32),
                        "first": np.array([[1.]] + [[0.] for _ in atoms[1:]], dtype=np.float32)}

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
        self._masses_are_set = False

        with CFG(filename) as cfg:
            self._manual_charges = {}
            for mol in cfg:
                self._mappings[mol.name] = []
                self._manual_charges[mol.name] = False
                molmap = self._mappings[mol.name]
                for name, typ, first, *atoms in mol:
                    charge = 0
                    try:
                        # Allow optional charge in mapping file
                        charge = float(first)
                        self._manual_charges[mol.name] = True
                    except ValueError:
                        atoms.insert(0, first)
                    assert atoms, "Bead {0} specification contains no atoms".format(name)
                    newbead = BeadMap(name=name, type=typ, atoms=atoms, charge=charge)
                    molmap.append(newbead)

        # TODO this only works with one moleculetype in one itp - extend this
        if itp is not None:
            with CFG(itp) as itp:
                atoms = {}
                for toks in itp["atoms"]:
                    # Store charge and mass
                    atoms[toks[4]] = (float(toks[6]), float(toks[7]))

                molname = itp["moleculetype"][0][0]
                for bead in self._mappings[molname]:
                    mass_array = np.array([[atoms[atom][1]] for atom in bead], dtype=np.float32)
                    bead.mass = sum(mass_array)
                    mass_array /= bead.mass
                    bead.weights["mass"] = mass_array

                    for atom in bead:
                        if self._manual_charges[molname]:
                            logger.warning("Charges assigned in mapping for molecule {0}, ignoring itp charges.".format(molname))
                        else:
                            bead.charge += atoms[atom][0]

                self._masses_are_set = True

    def __len__(self):
        return len(self._mappings)

    def __contains__(self, item):
        return item in self._mappings

    def __getitem__(self, item):
        return self._mappings[item]

    def __iter__(self):
        return iter(self._mappings)

    def _guess_atom_masses(self):
        """
        Guess atom masses from names
        """
        dist_dat_dir = os.path.join(dir_up(os.path.realpath(__file__), 2), "data")
        mass_file = os.path.join(dist_dat_dir, "atom_masses.json")
        with open(mass_file) as f:
            mass_dict = json.load(f)

        for mol_mapping in self._mappings.values():
            for bead in mol_mapping:
                if bead.mass == 0:
                    mass_array = np.zeros((len(bead.atoms), 1), dtype=np.float32)
                    for i, atom in enumerate(bead.atoms):
                        try:
                            mass = mass_dict[atom[:2]]
                        except KeyError:
                            try:
                                mass = mass_dict[atom[0]]
                            except KeyError:
                                msg = "Mass of atom {0} could not be automatically assigned, map_center=mass is not available."
                                raise RuntimeError(msg.format(atom))
                        mass_array[i] = mass
                    bead.mass = sum(mass_array)

                    if not np.all(mass_array):
                        msg = "Some atom masses could not be automatically assigned, map_center=mass is not available."
                        raise RuntimeError(msg)

                    mass_array /= bead.mass
                    bead.weights["mass"] = mass_array

        self._masses_are_set = True

    def _cg_frame_setup(self, aa_residues, name=None):
        """
        Create a new CG Frame and populate beads
        :param aa_residues: Iterable of atomistic residues to map from
        :param name: Name of Frame
        :return: New CG Frame instance
        """
        cgframe = Frame()
        cgframe.name = name

        missing_mappings = set()
        cg_bead_num = 1

        for aares in aa_residues:
            try:
                molmap = self._mappings[aares.name]
            except KeyError:
                if aares.name not in missing_mappings:
                    missing_mappings.add(aares.name)
                    logger.warning("A mapping has not been provided for '{0}' residues, they will not be mapped.".format(aares.name))
                continue

            cgres = Residue(name=aares.name, num=aares.num)
            cgres.atoms = [Atom(name=bmap.name, type=bmap.type, charge=bmap.charge, mass=bmap.mass, coords=np.zeros(3)) for bmap in molmap]

            for i, (bead, bmap) in enumerate(zip(cgres, molmap)):
                cgres.name_to_num[bead.name] = i
                bead.charge = bmap.charge
                bead.mass = bmap.mass
                bead.num = cg_bead_num

                cg_bead_num += 1

            cgframe.add_residue(cgres)
            cgframe.natoms += len(cgres)

        return cgframe

    def apply(self, frame, cgframe=None, exclude=None):
        """
        Apply the AA->CG mapping to an atomistic Frame.

        :param frame: Frame to which mapping will be applied
        :param cgframe: CG Frame to remap - optional
        :param exclude: Set of molecule names to exclude from mapping - e.g. solvent
        :return: Frame instance containing the CG frame
        """
        if self._map_center == "mass" and not self._masses_are_set:
            self._guess_atom_masses()

        if cgframe is None:
            # Frame needs initialising
            cgframe = self._cg_frame_setup(frame.residues, frame.name)

        cgframe.time = frame.time
        cgframe.number = frame.number
        cgframe.box = frame.box

        select_predicate = lambda res: res.name in self._mappings and not (exclude is not None and res.name in exclude)
        aa_residues = filter(select_predicate, frame)

        for aares, cgres in zip(aa_residues, cgframe):
            molmap = self._mappings[aares.name]

            for i, (bead, bmap) in enumerate(zip(cgres, molmap)):
                ref_coords = aares[bmap[0]].coords
                coords = np.array([aares[atom].coords for atom in bmap], dtype=np.float32)

                try:
                    weights = bmap.weights[self._map_center]
                except KeyError as e:
                    e.args = ("Unknown mapping type '{0}'.".format(e.args[0]),)
                    raise

                bead.coords = calc_coords_weight(ref_coords, coords, cgframe.box, weights)

        return cgframe


@numba.jit
def calc_coords_weight(ref_coords, coords, box, weights):
    """
    Calculate the coordinates of a single CG bead from weighted component atom coordinates.

    :param ref_coords: Coordinates of reference atom, usually first atom in bead
    :param coords: Array of coordinates of component atoms
    :param box: PBC box vectors
    :param weights: Array of atom weights, must sum to 1
    :return: Coordinates of CG bead
    """
    n = len(coords)
    result = np.zeros(3, dtype=np.float32)
    for i in range(n):
        tmp_coords = coords[i]
        tmp_coords -= ref_coords
        if box[0] * box[1] * box[2] != 0:
            tmp_coords -= box * np.rint(tmp_coords / box)
        result += weights[i] * tmp_coords
    result += ref_coords
    return result
