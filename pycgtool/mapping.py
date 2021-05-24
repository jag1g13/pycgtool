"""
This module contains classes required to perform an atomistic to coarse-grain mapping.

The Mapping class contains a dictionary of lists of BeadMaps.
Each list corresponds to a single molecule.
"""

import copy
import logging
import pathlib
import typing

import mdtraj
from mdtraj.formats.pdb import PDBTrajectoryFile
import numpy as np

from .frame import Frame
from .parsers import CFG, ITP

logger = logging.getLogger(__name__)


PathLike = typing.Union[str, pathlib.Path]


class BeadMap:
    """Class holding values relating to the AA->CG transformation for a single bead."""
    def __init__(self,
                 name: str,
                 num: int,
                 type: str = None,
                 atoms=None,
                 charge: float = 0,
                 mass: float = 0):
        """Create a single bead mapping.

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
        self.weights_dict = {
            "geom": np.array([1. / len(atoms) for _ in atoms],
                             dtype=np.float32),
            "first": np.array([1.] + [0. for _ in atoms[1:]],
                              dtype=np.float32)
        }
        self.weights = self.weights_dict["geom"]

    def __repr__(self):
        return f"BeadMap #{self.num} {self.name} type: {self.type} mass: {self.mass} charge: {self.charge}"

    def __iter__(self):
        """Iterate through the atom names from which the bead is made up.

        :return: Iterator over atoms
        """
        return iter(self.atoms)

    def __len__(self):
        return len(self.atoms)

    def __getitem__(self, item):
        return self.atoms[item]

    def guess_atom_masses(self) -> None:
        """Guess masses for the atoms inside this bead."""
        mass_array = np.zeros(len(self.atoms), dtype=np.float32)

        for i, atom in enumerate(self.atoms):
            try:
                mass = mdtraj.element.Element.getBySymbol(atom[:2]).mass

            except KeyError:
                try:
                    mass = mdtraj.element.Element.getBySymbol(atom[0]).mass

                except KeyError:
                    raise RuntimeError(
                        f"Mass of atom {atom} could not be automatically assigned, "
                        "map_center=mass is not available.")

            mass_array[i] = mass

        self.mass = sum(mass_array)

        if not np.all(mass_array):
            raise RuntimeError(
                "Some atom masses could not be automatically assigned, "
                "map_center=mass is not available")

        mass_array /= self.mass
        self.weights_dict["mass"] = mass_array


class VirtualMap(BeadMap):
    def __init__(self, name, num, type=None, atoms=None, charge=0):
        """Create a single bead mapping.

        :param str name: The name of the bead
        :param int num: The number of the bead
        :param str type: The bead type
        :param List[str] atoms: The CG bead names from which the bead position is determined
        :param float charge: The net charge on the bead
        """
        super().__init__(name,
                         num,
                         type=type,
                         atoms=atoms,
                         charge=charge,
                         mass=0.)

        self.gromacs_type_id_dict = {"geom": 1, "mass": 2}
        self.gromacs_type_id = self.gromacs_type_id_dict["geom"]

    def guess_atom_masses(self) -> None:
        """Virtual beads should not define atom masses."""
        return


class Mapping:
    """Class used to perform the AA->CG mapping.

    Contains a dictionary of lists of BeadMaps.  Each list corresponds to a single molecule.
    """
    def __init__(self, filename: PathLike, options, itp_filename: typing.Optional[PathLike] = None):
        """Read in the AA->CG mapping from a file.

        :param filename: File from which to read mapping
        :param options: Options from command line
        :param itp_filename: Optional GROMACS itp file for atom mass / charge
        """
        self._manual_charges = {}
        self._mappings = {}
        self._map_center = options.map_center
        self._virtual_map_center = options.virtual_map_center
        self._masses_are_set = False

        with CFG(filename) as cfg:
            for mol_name, mol_section in cfg.items():
                mol_map, manual_charges = self._mol_map_from_section(
                    mol_section)

                self._mappings[mol_name] = mol_map
                self._manual_charges[mol_name] = manual_charges

        if itp_filename is not None:
            self._load_itp(itp_filename)

        if ((self._map_center == "mass" or self._virtual_map_center == "mass")
                and not self._masses_are_set):
            self._guess_atom_masses()

        self._set_bead_weights()
        self._rename_atoms()

    def _load_itp(self, itp_filename: str) -> None:
        """Load mass and charge of atoms in mapping from a GROMACS itp topology file.

        :param itp_filename: ITP file to read
        """
        with ITP(itp_filename) as itp:
            for mol_name, mol_map in self._mappings.items():
                try:
                    itp_mol_entry = itp[mol_name]
                    manual_charges = self._manual_charges[mol_name]
                    if manual_charges:
                        logger.warning(
                            'Charges assigned in mapping for molecule %s, ignoring itp charges',
                            mol_name)

                    self._itp_section_into_mol_map(mol_map, itp_mol_entry,
                                                   manual_charges)

                except KeyError:
                    logger.warning(
                        "No itp information for molecule %s found in %s",
                        mol_name, itp.filename)

            self._masses_are_set = True

    def _set_bead_weights(self) -> None:
        """Set bead weights for the mapping center being used."""
        for mol_map in self._mappings.values():
            for bmap in mol_map:
                if isinstance(bmap, VirtualMap):
                    bmap.weights = bmap.weights_dict[self._virtual_map_center]
                    bmap.gromacs_type_id = bmap.gromacs_type_id_dict[
                        self._virtual_map_center]

                else:
                    bmap.weights = bmap.weights_dict[self._map_center]

    @staticmethod
    def _itp_section_into_mol_map(mol_map: typing.List[BeadMap], itp_mol_entry,
                                  manual_charges: bool) -> None:
        atoms = {}
        for toks in itp_mol_entry["atoms"]:
            # Store charge and mass
            atoms[toks[4]] = (float(toks[6]), float(toks[7]))

        for bead in mol_map:
            if not isinstance(bead, VirtualMap):
                mass_array = np.array([atoms[atom][1] for atom in bead])
                bead.mass = sum(mass_array)
                mass_array /= bead.mass
                bead.weights_dict["mass"] = mass_array

                if not manual_charges:
                    for atom in bead:
                        bead.charge += atoms[atom][0]

        for bead in mol_map:
            if isinstance(bead, VirtualMap):
                mass_array = np.array([
                    real_bead.mass
                    for real_bead in mol_map if real_bead.name in bead
                ])
                weights_array = mass_array / sum(mass_array)
                bead.weights_dict["mass"] = weights_array

                if not manual_charges:
                    bead.charge = sum([
                        real_bead.charge for real_bead in mol_map
                        if real_bead.name in bead
                    ])

    @staticmethod
    def _mol_map_from_section(
            mol_section) -> typing.Tuple[typing.List[BeadMap], bool]:
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
                    raise ValueError(f'Invalid line prefix "{name}" in mapping file') from exc

            try:
                # Allow optional charge in mapping file
                # Using this even once in a molecule enables it for the whole molecule
                charge = float(first)
                manual_charges = True

            except ValueError:
                # Is an atom name, not a charge
                atoms.insert(0, first)

            if not atoms:
                raise ValueError(f'Bead {name} specification contains no atoms')

            mol_map.append(
                bead_class(name, i, type=typ, atoms=atoms, charge=charge))

        return mol_map, manual_charges

    def _rename_atoms(self) -> None:
        """Rename residues and atoms in mapping according to MDTraj conventions.

        This means that users can create mappings using the same names as are
        in the input topology file and not have this broken by MDTraj renaming
        the residues and atoms.
        """
        # We need to use a couple of protected methods from MDTraj here
        # pylint: disable=protected-access

        PDBTrajectoryFile._loadNameReplacementTables()

        atom_replacements = {}
        for mol_name in list(self._mappings.keys()):
            mol_map = copy.deepcopy(self._mappings[mol_name])

            try:
                mol_name = PDBTrajectoryFile._residueNameReplacements[mol_name]
                atom_replacements = PDBTrajectoryFile._atomNameReplacements[
                    mol_name]

            except KeyError:
                continue

            self._mappings[mol_name] = mol_map

            for bead in mol_map:
                bead.atoms = [
                    atom_replacements[atom]
                    if atom in atom_replacements else atom
                    for atom in bead.atoms
                ]

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

    def _guess_atom_masses(self) -> None:
        """Guess atom masses from their names."""
        for mol_mapping in self._mappings.values():
            for bead in mol_mapping:
                if not bead.mass:
                    bead.guess_atom_masses()

            # Set virtual bead masses
            # Do this afterwards as it depends on real atom masses
            for bead in mol_mapping:
                if isinstance(bead, VirtualMap):
                    mass_array = np.array([
                        real_bead.mass
                        for real_bead in mol_mapping if real_bead.name in bead
                    ], dtype=np.float32)  # yapf: disable

                    weights_array = mass_array / sum(mass_array)
                    bead.weights_dict["mass"] = weights_array

        self._masses_are_set = True

    def _cg_frame_setup(
            self, aa_residues: typing.Iterable[mdtraj.core.topology.Residue]):
        """Create a new CG Frame and populate beads
        :param aa_residues: Iterable of atomistic residues to map from
        :return: New CG Frame instance
        """
        logger.info('Initialising output frame')
        cg_frame = Frame()
        missing_mappings = set()

        for aa_res in aa_residues:
            try:
                mol_map = self._mappings[aa_res.name]

            except KeyError:
                if aa_res.name not in missing_mappings:
                    missing_mappings.add(aa_res.name)
                    logger.warning(
                        "A mapping has not been provided for '%s' residues, they will not be mapped",
                        aa_res.name)

                continue

            cg_res = cg_frame.add_residue(aa_res.name)
            for bmap in mol_map:
                cg_frame.add_atom(bmap.name, None, cg_res)

        logger.info('Finished initialising output frame')
        return cg_frame

    def apply(self, frame: Frame, cg_frame: typing.Optional[Frame] = None):
        """Apply the AA->CG mapping to an atomistic Frame.

        :param frame: Frame to which mapping will be applied
        :param cgframe: CG Frame to remap - optional
        :return: Frame instance containing the CG frame
        """
        if cg_frame is None:
            # Frame needs initialising
            cg_frame = self._cg_frame_setup(frame.residues)

        cg_frame.time = frame.time
        cg_frame.unitcell_lengths = frame.unitcell_lengths
        cg_frame.unitcell_angles = frame.unitcell_angles

        unitcell_lengths = cg_frame.unitcell_lengths
        if not np.all(unitcell_lengths):
            unitcell_lengths = None

        logger.info('Applying AA->CG mapping')
        residues_to_map = (res for res in frame.residues
                           if res.name in self._mappings)
        for aa_res, cg_res in zip(residues_to_map, cg_frame.residues):
            mol_map = self._mappings[aa_res.name]

            for bead, bmap in zip(cg_res.atoms, mol_map):
                if isinstance(bmap, VirtualMap):
                    continue

                ref_coords = aa_res.atom(bmap[0]).coords
                if len(bmap) == 1:
                    bead.coords = ref_coords

                else:
                    coords = np.array(
                        [aa_res.atom(atom).coords for atom in bmap],
                        dtype=np.float32)
                    bead.coords = calc_coords_weight(ref_coords, coords,
                                                     bmap.weights,
                                                     unitcell_lengths)

            for bead, bmap in zip(cg_res.atoms, mol_map):
                if isinstance(bmap, VirtualMap):
                    coords = np.array(
                        [cg_res.atom(atom).coords for atom in bmap],
                        dtype=np.float32)
                    bead.coords = calc_coords_weight(ref_coords, coords,
                                                     bmap.weights,
                                                     unitcell_lengths)

        cg_frame.build_trajectory()

        logger.info('Finished applying AA->CG mapping')
        return cg_frame


def calc_coords_weight(ref_coords, coords, weights, box=None):
    """Calculate the coordinates of a single CG bead from weighted component atom coordinates.

    :param ref_coords: Coordinates of reference atom, usually first atom in bead
    :param coords: Array of coordinates of component atoms
    :param weights: Array of atom weights, must sum to 1
    :param box: PBC box vectors
    :return: Coordinates of CG bead
    """
    vectors = coords - ref_coords
    if box is not None:
        vectors -= box * np.rint(vectors / box)

    # Reshape weights array to match atom positions
    weights_t = weights[np.newaxis, np.newaxis].T
    result = np.sum(weights_t * vectors, axis=0)
    result += ref_coords
    return result
