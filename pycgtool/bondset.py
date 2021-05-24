"""Module containing classes to calculate bonded properties from a Frame.

BondSet contains a dictionary of lists of Bonds.  Each list corresponds to a single molecule.
"""

import itertools
import logging
import math
import pathlib
import typing

import numpy as np
import mdtraj

from .mapping import VirtualMap
from .functionalforms import get_functional_forms
from .parsers.cfg import CFG
from .util import (
    circular_mean,
    circular_variance,
    extend_graph_chain,
    file_write_lines,
    sliding,
    transpose_and_sample,
)

logger = logging.getLogger(__name__)

PathLike = typing.Union[str, pathlib.Path]


class Bond:
    """Class holding the properties of a single bonded term.

    Bond lengths, angles and dihedrals are all equivalent, distinguished by the number of atoms present.
    """
    def __init__(self,
                 atoms: typing.Iterable[str],
                 atom_numbers: typing.Optional[typing.Iterable[int]] = None,
                 func_form=None):
        """Create a single bond definition.

        :param List[str] atoms: List of atom names defining the bond
        :param List[int] atom_numbers: List of atom numbers defining the bond
        :param func_form: Functional form to use for Boltzmann Inversion
        """
        self.atoms = atoms
        self.atom_numbers = atom_numbers
        self.values = []
        self.eqm = None
        self.fconst = None

        self._func_form = func_form
        self.gromacs_type_id = func_form.gromacs_type_id_by_natoms(len(atoms))

    def __iter__(self):
        return iter(self.atoms)

    def boltzmann_invert(self, temp: float = 310) -> None:
        """Perform Boltzmann Inversion using measured values of bond to calculate equilibrium value and force constant.

        :param temp: Temperature at which the simulation was performed
        """
        if len(self.values) == 0:
            raise ValueError("No bonds were measured between atoms {0}".format(
                self.atoms))

        values = np.array(self.values)

        with np.errstate(divide="raise"):
            self.eqm = self._func_form.eqm(values, temp)

            # TODO: consider moving this correction to the functional form
            # If dihedral, get value to shift cosine function by, NOT equilibrium value
            if len(self.atoms) == 4:
                two_pi = 2 * np.pi
                self.eqm -= np.pi
                self.eqm = (((self.eqm + np.pi) % two_pi) - np.pi)

            try:
                self.fconst = self._func_form.fconst(values, temp)

            except FloatingPointError:
                # Happens when variance is 0, i.e. we only have one value
                self.fconst = float("inf")

    def __repr__(self) -> str:
        try:
            return "<Bond containing atoms {0} with r_0 {1:.3f} and force constant {2:.3e}>".format(
                ", ".join(self.atoms), self.eqm, self.fconst)
        except (AttributeError, TypeError):
            return "<Bond containing atoms {0}>".format(", ".join(self.atoms))


class BondSet:
    """Class used to perform bond measurements in a Frame.

    BondSet contains a dictionary of lists of Bonds. Each list corresponds to a single molecule.
    """
    def _get_default_func_forms(self, options):

        try:
            default_fc = options.default_fc

        except AttributeError:
            default_fc = False

        if default_fc:
            default_forms = [
                'MartiniDefaultLength', 'MartiniDefaultAngle',
                'MartiniDefaultDihedral'
            ]
        else:
            default_forms = ['Harmonic', 'CosHarmonic', 'Harmonic']

        functional_forms_map = get_functional_forms()

        functional_forms = [None, None]
        functional_forms.extend(
            [functional_forms_map[ff].value for ff in default_forms])

        try:
            functional_forms[2] = functional_forms_map[
                options.length_form].value
        except AttributeError:
            pass

        try:
            functional_forms[3] = functional_forms_map[
                options.angle_form].value
        except AttributeError:
            pass

        try:
            functional_forms[4] = functional_forms_map[
                options.dihedral_form].value
        except AttributeError:
            pass

        return functional_forms

    def __init__(self, filename: PathLike, options):
        """Read in bonds from a file.

        :param filename: File to read
        :return: Instance of BondSet
        """
        self._molecules = {}

        self._fconst_constr_threshold = options.constr_threshold

        try:
            self._temperature = options.temperature
        except AttributeError:
            self._temperature = 310.

        self._functional_forms = self._get_default_func_forms(options)

        with CFG(filename) as cfg:
            for mol_name, mol_section in cfg.items():
                self._molecules[mol_name] = []
                mol_bonds = self._molecules[mol_name]

                angles_defined = False
                for atomlist in mol_section:
                    is_angle_or_dihedral = len(atomlist) > 2
                    mean_function = circular_mean if is_angle_or_dihedral else np.nanmean
                    variance_function = circular_variance if is_angle_or_dihedral else np.nanvar

                    # Construct instance of Boltzmann Inversion function and
                    # inject dependencies for mean and variance functions
                    # TODO: Should we allow overriding functional forms per bond?
                    func_form = self._functional_forms[len(atomlist)](
                        mean_function, variance_function)

                    if {x for x in atomlist if atomlist.count(x) > 1}:
                        raise ValueError(
                            "Defined bond '{0}' contains duplicate atoms".
                            format(atomlist))

                    mol_bonds.append(Bond(atoms=atomlist, func_form=func_form))
                    if len(atomlist) > 2:
                        angles_defined = True

                if not angles_defined:
                    angles, dihedrals = self._create_angles(mol_bonds)

                    if options.generate_angles:
                        for atomlist in angles:
                            mol_bonds.append(
                                Bond(atoms=atomlist,
                                     func_form=self._functional_forms[3](
                                         circular_mean, circular_variance)))

                    if options.generate_dihedrals:
                        for atomlist in dihedrals:
                            mol_bonds.append(
                                Bond(atoms=atomlist,
                                     func_form=self._functional_forms[4](
                                         circular_mean, circular_variance)))

    @staticmethod
    def _create_angles(mol_bonds):
        """Create angles and dihedrals from bonded topology.

        :param mol_bonds: List of bonds within a molecule to generate angles for
        :return: List of angle and dihedral bond name tuples
        """
        bonds = [bond.atoms for bond in mol_bonds if len(bond.atoms) == 2]
        angles = extend_graph_chain(bonds, bonds)
        dihedrals = extend_graph_chain(angles, bonds)
        return angles, dihedrals

    def get_bonds(
        self,
        mol: str,
        natoms: int,
        select: typing.Callable[[Bond], bool] = lambda x: True
    ) -> typing.List[Bond]:
        """Return list of bonds from molecule containing natoms atoms.

        :param str mol: Molecule name
        :param int natoms: Number of atoms in bond, i.e. 2 for bond length, 3 for angle, 4 for dihedral
        :param function select: Optional lambda, return only bonds for which this is True
        :return List[Bond]: List of bonds
        """
        if natoms == -1:
            return [bond for bond in self._molecules[mol] if select(bond)]

        return [
            bond for bond in self._molecules[mol]
            if len(bond.atoms) == natoms and select(bond)
        ]

    def get_bond_lengths(self, mol, with_constr=False):
        """Return list of all bond lengths in molecule.  May include constraints.

        :param mol: Molecule name to return bonds for
        :param with_constr: Include constraints?
        :return: List of bonds
        """
        if with_constr:
            return [
                bond for bond in self._molecules[mol] if len(bond.atoms) == 2
            ]

        return [
            bond for bond in self._molecules[mol] if len(bond.atoms) == 2
            and bond.fconst < self._fconst_constr_threshold
        ]

    def get_bond_length_constraints(self, mol):
        """Return list of all bond length constraints in molecule.

        :param mol: Molecule name to return bonds for
        :return: List of bonds
        """
        return [
            bond for bond in self._molecules[mol] if len(bond.atoms) == 2
            and bond.fconst >= self._fconst_constr_threshold
        ]

    def get_bond_angles(self, mol, exclude_triangle=True):
        """Return list of all bond angles in molecule.

        :param mol: Molecule name to return bonds for
        :param exclude_triangle: Exclude angles that are part of a triangle?
        :return: List of bonds
        """
        angles = [
            bond for bond in self._molecules[mol] if len(bond.atoms) == 3
        ]

        if exclude_triangle:
            edges = [
                tuple(bond.atoms)
                for bond in self.get_bond_lengths(mol, with_constr=True)
            ]

            def is_triangle(atoms):
                triangle_edges = 0
                for j in range(3):
                    if (atoms[j - 1],
                            atoms[j]) in edges or (atoms[j],
                                                   atoms[j - 1]) in edges:
                        triangle_edges += 1
                return triangle_edges >= 3

            angles = [
                angle for angle in angles if not is_triangle(angle.atoms)
            ]

        return angles

    def get_bond_dihedrals(self, mol):
        """Return list of all bond dihedrals in molecule.

        :param mol: Molecule name to return bonds for
        :return: List of bonds
        """
        return [bond for bond in self._molecules[mol] if len(bond.atoms) == 4]

    def get_virtual_beads(self, mapping):
        """Return list of all virtual beads in molecule
        :param mapping:
        :return: list of virtual beads
        """
        return [bead for bead in mapping if isinstance(bead, VirtualMap)]

    def _populate_atom_numbers(self, mapping):
        """Add atom numbers to all bonds.

        Uses previously defined atom names.

        :param mapping: AA->CG mapping used to collect atom/bead numbers
        """
        for mol in self._molecules:
            try:
                molmap = mapping[mol]
            except KeyError:
                # Molecule was not mapped, can't create itp
                continue

            index = [bead.name for bead in molmap]
            for bond in self._molecules[mol]:
                # TODO this causes issue #8
                try:
                    bond.atom_numbers = [
                        index.index(atom.lstrip("+-")) for atom in bond.atoms
                    ]
                except ValueError as e:
                    missing = [
                        atom for atom in bond.atoms
                        if atom.lstrip("+-") not in index
                    ]
                    e.args = (
                        "Bead(s) {0} do(es) not exist in residue {1}".format(
                            missing, mol), )
                    raise

    def write_itp(self, filename, mapping):
        """Output a GROMACS .itp file containing atoms/beads and bonded terms.

        :param filename: Name of output file
        :param mapping: AA->CG Mapping from which to collect bead properties
        """
        self._populate_atom_numbers(mapping)
        file_write_lines(filename, self.itp_text(mapping))

    def itp_text(self, mapping):
        atom_template = {
            "nomass": "{0:4d} {1:4s} {2:4d} {3:4s} {4:4s} {5:4d} {6:8.3f}",
            "mass":
            "{0:4d} {1:4s} {2:4d} {3:4s} {4:4s} {5:4d} {6:8.3f} {7:8.3f}"
        }

        def write_bond_angle_dih(bonds,
                                 section_header,
                                 print_fconst=True,
                                 multiplicity=None,
                                 rad2deg=False):
            ret_lines = []
            if bonds:
                ret_lines.append("\n[ {0:s} ]".format(section_header))
            for bond in bonds:
                line = " ".join([
                    "{0:4d}".format(atnum + 1) for atnum in bond.atom_numbers
                ])
                eqm = math.degrees(bond.eqm) if rad2deg else bond.eqm
                line += " {0:4d} {1:12.5f}".format(bond.gromacs_type_id, eqm)
                if print_fconst:
                    line += " {0:12.5f}".format(bond.fconst)
                if multiplicity is not None:
                    line += " {0:4d}".format(multiplicity)
                ret_lines.append(line)

            return ret_lines

        ret_lines = [
            ";",
            "; Topology prepared automatically using PyCGTOOL",
            "; James Graham <J.A.Graham@soton.ac.uk> 2016",
            "; University of Southampton",
            "; https://github.com/jag1g13/pycgtool",
            ";"
        ]  # yapf: disable

        # Print molecule
        for mol in self._molecules:
            if mol not in mapping:
                logger.warning((
                    "Molecule '%s' present in bonding file, but not in mapping. "
                    "Parameters have not been calculated."), mol)
                continue

            if not all(bond.fconst is not None for bond in self._molecules[mol]):
                logger.warning(("Molecule '%s' has no measured bond values. "
                                "Parameters have not been calculated."), mol)
                continue

            ret_lines.append("\n[ moleculetype ]")
            ret_lines.append("{0:4s} {1:4d}".format(mol, 1))

            ret_lines.append("\n[ atoms ]")

            for i, bead in enumerate(mapping[mol], start=1):
                # atnum type resnum resname atname c-group charge (mass)
                if isinstance(bead, VirtualMap):
                    ret_lines.append(atom_template["mass"].format(
                        i, bead.type, 1, mol, bead.name, i, bead.charge,
                        bead.mass))
                else:
                    ret_lines.append(atom_template["nomass"].format(
                        i, bead.type, 1, mol, bead.name, i, bead.charge))

            virtual_beads = self.get_virtual_beads(mapping[mol])
            if len(virtual_beads) != 0:
                ret_lines.append("\n[ virtual_sitesn ]")
                excl_lines = ["\n[ exclusions ]"
                              ]  # Exclusions section for virtual sites

                for vbead in virtual_beads:
                    cg_ids = sorted([
                        bead.num + 1 for bead in mapping[mol]
                        if bead.name in vbead.atoms
                    ])
                    cg_ids_string = " ".join(map(str, cg_ids))
                    ret_lines.append("{0:^6d} {1:^6d} {2}".format(
                        vbead.num + 1, vbead.gromacs_type_id, cg_ids_string))
                    vsite_exclusions = "{} ".format(vbead.num +
                                                    1) + cg_ids_string
                    excl_lines.append(vsite_exclusions)

                ret_lines.extend(excl_lines)

            ret_lines.extend(
                write_bond_angle_dih(self.get_bond_lengths(mol), "bonds"))
            ret_lines.extend(
                write_bond_angle_dih(self.get_bond_angles(mol),
                                     "angles",
                                     rad2deg=True))
            ret_lines.extend(
                write_bond_angle_dih(self.get_bond_dihedrals(mol),
                                     "dihedrals",
                                     multiplicity=1,
                                     rad2deg=True))
            ret_lines.extend(
                write_bond_angle_dih(self.get_bond_length_constraints(mol),
                                     "constraints",
                                     print_fconst=False))

        return ret_lines

    def apply(self, frame):
        """Calculate bond lengths/angles for all trajectory frames with a Frame instance and store into Bonds.

        :param frame: Frame from which to calculate values
        """
        calc = {
            2: mdtraj.compute_distances,
            3: mdtraj.compute_angles,
            4: mdtraj.compute_dihedrals
        }

        for mol_name, mol_bonds in self._molecules.items():
            for bond in mol_bonds:
                bond_indices = []

                for prev_residue, residue, next_residue in sliding(
                        frame.residues):
                    adj_res = {"-": prev_residue, "+": next_residue}

                    # Need access to adjacent residues, so can't just iterate over these directly
                    if residue.name == mol_name:
                        try:
                            bond_indices.append([
                                adj_res.get(atom_name[0], residue).atom(
                                    atom_name.lstrip('-+')).index
                                for atom_name in bond.atoms
                            ])

                        except AttributeError:
                            # AttributeError is raised when residues on end of chain calc bond to next
                            pass

                vals = calc[len(bond.atoms)](frame._trajectory, bond_indices)
                bond.values = vals.flatten()

    def boltzmann_invert(self) -> None:
        """Perform Boltzmann Inversion of all bonds to calculate equilibrium value and force constant."""
        for bond in itertools.chain(*self._molecules.values()):
            try:
                bond.boltzmann_invert(temp=self._temperature)
            except ValueError:
                pass

    @staticmethod
    def _get_lines_for_bond_dump(bonds,
                                 target_number: typing.Optional[int] = None,
                                 rad2deg: bool = False) -> typing.List[str]:
        """Return a dump of bond measurements as a list of lines of text.

        :param bonds: Iterable of bonds
        :param int target_number: Sample size of bonds to dump - None will dump all bonds
        :param bool rad2deg: Should bond measurements be converted from radians to degrees?
        :return List[str]: Lines of bond measurements dump
        """
        ret_lines = []
        for row in transpose_and_sample((bond.values for bond in bonds),
                                        n=target_number):
            if rad2deg:
                row = [math.degrees(val) for val in row]

            ret_lines.append((len(row) * '{:12.5f}').format(*row))
        return ret_lines

    def dump_values(self, target_number: int = 10000,
                    output_dir: typing.Optional[typing.Union[pathlib.Path, str]] = None) -> None:
        """Output measured bond values to files for length, angles and dihedrals.

        :param int target_number: Approx number of sample measurements to output. If None, all samples will be output.
        :param output_dir: Directory to write output files to.
        """
        if output_dir is None:
            output_dir = '.'

        # Cast to Path type
        output_dir = pathlib.Path(output_dir)

        for mol in self._molecules:
            if mol == 'SOL':
                continue

            bonds = self.get_bond_lengths(mol, with_constr=True)
            if bonds:
                lines = BondSet._get_lines_for_bond_dump(bonds, target_number)
                file_write_lines(output_dir.joinpath(f'{mol}_length.dat'), lines)

            bonds = self.get_bond_angles(mol)
            if bonds:
                lines = BondSet._get_lines_for_bond_dump(bonds,
                                                         target_number,
                                                         rad2deg=True)
                file_write_lines(output_dir.joinpath(f'{mol}_angle.dat'), lines)

            bonds = self.get_bond_dihedrals(mol)
            if bonds:
                lines = BondSet._get_lines_for_bond_dump(bonds,
                                                         target_number,
                                                         rad2deg=True)
                file_write_lines(output_dir.joinpath(f'{mol}_dihedral.dat'), lines)

    def __len__(self):
        return len(self._molecules)

    def __contains__(self, item):
        return item in self._molecules

    def __getitem__(self, item):
        return self._molecules[item]

    def __iter__(self):
        return iter(self._molecules)
