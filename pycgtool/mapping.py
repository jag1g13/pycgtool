"""
This module contains classes required to perform an atomistic to coarse-grain mapping.

The Mapping class contains a dictionary of lists of BeadMaps.
Each list corresponds to a single molecule.
"""

import numpy as np
import logging
import json
import os

from .frame import Frame
from .parsers import CFG, ITP
from .util import dir_up

try:
    import numba
except ImportError:
    from .util import NumbaDummy
    numba = NumbaDummy()

logger = logging.getLogger(__name__)


class BeadMap:
    """
    POD class holding values relating to the AA->CG transformation for a single bead.
    """
    __slots__ = ["name", "num", "type", "atoms", "charge", "mass", "weights", "weights_dict"]

    def __init__(self, name, num, type=None, atoms=None, charge=0, mass=0):
        """
        Create a single bead mapping.

        :param str name: The name of the bead
        :param int num: The number of the bead
        :param str type: The bead type
        :param List[str] atoms: The atom names from which the bead is made up
        :param float charge: The net charge on the bead
        :param float mass: The total bead mass
        """
        self.name = name
        self.num = num
        self.type = type
        self.mass = mass
        self.charge = charge

        self.atoms = atoms
        # NB: Mass weights are added in Mapping.__init__ if an itp file is provided
        self.weights_dict = {"geom": np.array([[1. / len(atoms)] for _ in atoms], dtype=np.float32),
                             "first": np.array([[1.]] + [[0.] for _ in atoms[1:]], dtype=np.float32)}
        self.weights = self.weights_dict["geom"]

    def __repr__(self):
        return "BeadMap #{0} {1} type: {2} mass: {3} charge: {4}".format(
            self.num, self.name, self.type, self.mass, self.charge
        )

    def add_missing_data(self, other):
        assert self.name == other.name
        assert self.num == other.num

        for attr in ("type", "mass", "charge"):
            if getattr(self, attr) is None:
                setattr(self, attr, getattr(other, attr))

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


class VirtualMap(BeadMap):
    __slots__ = ["name", "type", "atoms", "charge", "mass", "weights", "weights_dict",
                 "gromacs_type_id_dict", "gromacs_type_id"]

    def __init__(self, name, num, type=None, atoms=None, charge=0):
        """
        Create a single bead mapping.

        :param str name: The name of the bead
        :param int num: The number of the bead
        :param str type: The bead type
        :param List[str] atoms: The CG bead names from which the bead position is determined
        :param float charge: The net charge on the bead
        """
        super().__init__(name, num, type=type, atoms=atoms, charge=charge, mass=0.)

        self.gromacs_type_id_dict = {"geom": 1, "mass": 2}
        self.gromacs_type_id = self.gromacs_type_id_dict["geom"]


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
        self._manual_charges = {}
        self._mappings = {}
        self._map_center = options.map_center
        self._virtual_map_center = options.virtual_map_center
        self._masses_are_set = False

        with CFG(filename) as cfg:
            for mol_name, mol_section in cfg.items():
                mol_map, manual_charges = self._mol_map_from_section(mol_section)

                self._mappings[mol_name] = mol_map
                self._manual_charges[mol_name] = manual_charges

        if itp is not None:
            with ITP(itp) as itp:
                for molname in self._mappings:
                    try:
                        molentry = itp[molname]
                        atoms = {}
                        for toks in molentry["atoms"]:
                            # Store charge and mass
                            atoms[toks[4]] = (float(toks[6]), float(toks[7]))
                        for bead in self._mappings[molname]:
                            if not isinstance(bead, VirtualMap):
                                mass_array = np.array([[atoms[atom][1]] for atom in bead], dtype=np.float32)
                                bead.mass = sum(mass_array)
                                print(mass_array, sum(mass_array))
                                mass_array /= bead.mass
                                bead.weights_dict["mass"] = mass_array

                                for atom in bead:
                                    if self._manual_charges[molname]:
                                        logger.warning("Charges assigned in mapping for molecule {0}, "
                                                       "ignoring itp charges.".format(molname))
                                    else:
                                        bead.charge += atoms[atom][0]

                        for bead in self._mappings[molname]:
                            if isinstance(bead, VirtualMap):
                                mass_array = np.array([real_bead.mass for real_bead in self._mappings[molname]
                                                       if real_bead.name in bead], dtype=np.float32)
                                weights_array = mass_array / sum(mass_array)
                                bead.weights_dict["mass"] = weights_array
                                if self._manual_charges[molname]:
                                    logger.warning("Charges assigned in mapping for molecule {0}, "
                                                   "ignoring itp charges.".format(molname))
                                else:
                                    charges = [real_bead.charge for real_bead in self._mappings[molname]
                                               if real_bead.name in bead]
                                    bead.charge = sum(charges)
                    except KeyError:
                        logger.warning("No itp information for molecule {0} found in {1}".format(molname, itp.filename))

                self._masses_are_set = True

        if ((self._map_center == "mass" or self._virtual_map_center == "mass")
                and not self._masses_are_set):
            self._guess_atom_masses()

        for molname, mapping in self._mappings.items():
            for bmap in mapping:
                if isinstance(bmap, VirtualMap):
                    bmap.weights = bmap.weights_dict[self._virtual_map_center]
                    bmap.gromacs_type_id = bmap.gromacs_type_id_dict[self._virtual_map_center]

                else:
                    bmap.weights = bmap.weights_dict[self._map_center]

    @staticmethod
    def _mol_map_from_section(mol_section):
        mol_map = []
        manual_charges = False

        prefix_to_class = {'@v': VirtualMap}

        for i, (name, typ, first, *atoms) in enumerate(mol_section):
            bead_class = BeadMap
            charge = 0

            if name.startswith('@'):
                try:
                    bead_class = prefix_to_class[name]
                    name, typ, first, *atoms = mol_section[i][1:]

                except KeyError as exc:
                    raise ValueError(
                        '"{0}" line prefix invalid'.format(name)) from exc

            try:
                # Allow optional charge in mapping file
                # Using this even once in a molecule enables it for the whole molecule
                charge = float(first)
                manual_charges = True

            except ValueError:
                # Is an atom name, not a charge
                atoms.insert(0, first)

            if not atoms:
                # TODO should this stop execution?
                logger.warning('Bead %s specification contains no atoms', name)

            mol_map.append(
                bead_class(name, i, type=typ, atoms=atoms, charge=charge))

        return mol_map, manual_charges

    def __len__(self):
        return len(self._mappings)

    def __contains__(self, item):
        return item in self._mappings

    def __getitem__(self, item):
        return self._mappings[item]

    def __iter__(self):
        return iter(self._mappings)

    def items(self):
        return self._mappings.items()

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
                if bead.mass == 0 and not isinstance(bead, VirtualMap):
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
                    bead.weights_dict["mass"] = mass_array

            #set virtual bead masses#
            for bead in mol_mapping:
                if isinstance(bead, VirtualMap):
                    mass_array = np.array([real_bead.mass for real_bead in mol_mapping if real_bead.name in bead], dtype=np.float32)
                    weights_array = mass_array / sum(mass_array)
                    bead.weights_dict["mass"] = weights_array

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
            cgres.atoms = [Atom(bmap.name, bmap.num, type=bmap.type, charge=bmap.charge, mass=bmap.mass, coords=np.zeros(3)) for bmap in molmap]

            for i, (bead, bmap) in enumerate(zip(cgres, molmap)):
                cgres.name_to_num[bead.name] = i
                bead.charge = bmap.charge
                bead.mass = bmap.mass
                bead.num = cg_bead_num

                cg_bead_num += 1

            cgframe.add_residue(cgres)
            cgframe.natoms += len(cgres)

        return cgframe

    def apply(self, frame, cgframe=None):
        """
        Apply the AA->CG mapping to an atomistic Frame.

        :param frame: Frame to which mapping will be applied
        :param cgframe: CG Frame to remap - optional
        :return: Frame instance containing the CG frame
        """
        if cgframe is None:
            # Frame needs initialising
            cgframe = self._cg_frame_setup(frame.yield_resname_in(self._mappings), frame.name)

        cgframe.time = frame.time
        cgframe.number = frame.number
        cgframe.box = frame.box

        coord_func = calc_coords_weight if frame.box[0] * frame.box[1] * frame.box[2] else calc_coords_weight_nobox

        for aares, cgres in zip(frame.yield_resname_in(self._mappings), cgframe):
            molmap = self._mappings[aares.name]

            virtual_beads= []
            virtual_bmap = []
            for i, (bead, bmap) in enumerate(zip(cgres, molmap)):
                if isinstance(bmap, VirtualMap):
                    virtual_beads.append(bead)
                    virtual_bmap.append(bmap)
                    continue

                ref_coords = aares[bmap[0]].coords
                if len(bmap) == 1:
                    bead.coords = ref_coords
                    continue

                coords = np.asarray([aares[atom].coords for atom in bmap], dtype=np.float32)
                bead.coords = coord_func(ref_coords, coords, cgframe.box, bmap.weights)

            for bead, bmap in zip(virtual_beads, virtual_bmap):
                coords = np.asarray([cgres[atom].coords for atom in bmap], dtype=np.float32)
                bead.coords = coord_func(ref_coords, coords, cgframe.box, bmap.weights)

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
    vectors = coords - ref_coords
    vectors -= box * np.rint(vectors / box)
    result = np.sum(weights * vectors, axis=0)
    result += ref_coords
    return result


@numba.jit
def calc_coords_weight_nobox(ref_coords, coords, box, weights):
    """
    Calculate the coordinates of a single CG bead from weighted component atom coordinates.

    :param ref_coords: Coordinates of reference atom, usually first atom in bead
    :param coords: Array of coordinates of component atoms
    :param box: PBC box vectors
    :param weights: Array of atom weights, must sum to 1
    :return: Coordinates of CG bead
    """
    vectors = coords - ref_coords
    result = np.sum(weights * vectors, axis=0)
    result += ref_coords
    return result
