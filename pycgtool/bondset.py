"""
Module containing classes to calculate bonded properties from a Frame.

BondSet contains a dictionary of lists of Bonds.  Each list corresponds to a single molecule.
"""

import itertools
import math
import logging

import numpy as np


try:
    from tqdm import tqdm
except ImportError:
    from .util import tqdm_dummy as tqdm

from .util import sliding, dist_with_pbc, transpose_and_sample
from .util import extend_graph_chain, backup_file
from .util import vector_len, vector_cross, vector_angle, vector_angle_signed
from .parsers.cfg import CFG
from .functionalforms import FunctionalForms

logger = logging.getLogger(__name__)


class Bond:
    """
    Class holding the properties of a single bonded term.

    Bond lengths, angles and dihedrals are all equivalent, distinguished by the number of atoms present.
    """
    __slots__ = ["atoms", "atom_numbers", "values", "eqm", "fconst", "gromacs_type_id", "_func_form"]

    def __init__(self, atoms, atom_numbers=None, func_form=None):
        """
        Create a single bond definition.

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

    def __len__(self):
        return len(self.atoms)

    def __iter__(self):
        return iter(self.atoms)

    def boltzmann_invert(self, temp=310):
        """
        Perform Boltzmann Inversion using measured values of bond to calculate equilibrium value and force constant.

        :param temp: Temperature at which the simulation was performed
        """
        if not self.values:
            raise ValueError("No bonds were measured between atoms {0}".format(self.atoms))

        values = np.array(self.values)

        with np.errstate(divide="raise"):
            self.eqm = self._func_form.eqm(values, temp)
            try:
                self.fconst = self._func_form.fconst(values, temp)
            except FloatingPointError:
                # Happens when variance is 0, i.e. we only have one value
                self.fconst = float("inf")

    def __repr__(self):
        try:
            return "<Bond containing atoms {0} with r_0 {1:.3f} and force constant {2:.3e}>".format(
                ", ".join(self.atoms), self.eqm, self.fconst
            )
        except (AttributeError, TypeError):
            return "<Bond containing atoms {0}>".format(", ".join(self.atoms))


class BondSet:
    """
    Class used to perform bond measurements in a Frame.

    BondSet contains a dictionary of lists of Bonds.  Each list corresponds to a single molecule.
    """
    def __init__(self, filename, options):
        """
        Read in bonds from a file.

        :param filename: File to read
        :return: Instance of BondSet
        """
        self._molecules = {}

        self._fconst_constr_threshold = options.constr_threshold

        try:
            self._temperature = options.temperature
        except AttributeError:
            self._temperature = 310.

        try:
            self._default_fc = options.default_fc
        except AttributeError:
            self._default_fc = False

        # Setup default functional forms
        functional_forms = FunctionalForms()
        if self._default_fc:
            default_forms = ["MartiniDefaultLength", "MartiniDefaultAngle", "MartiniDefaultDihedral"]
        else:
            default_forms = ["Harmonic", "CosHarmonic", "Harmonic"]
        self._functional_forms = [None, None]
        self._functional_forms.extend(map(lambda x: functional_forms[x], default_forms))

        try:
            self._functional_forms[2] = functional_forms[options.length_form]
        except AttributeError:
            pass

        try:
            self._functional_forms[3] = functional_forms[options.angle_form]
        except AttributeError:
            pass

        try:
            self._functional_forms[4] = functional_forms[options.dihedral_form]
        except AttributeError:
            pass

        with CFG(filename) as cfg:
            for mol_name, mol_section in cfg.items():
                self._molecules[mol_name] = []
                mol_bonds = self._molecules[mol_name]

                angles_defined = False
                for atomlist in mol_section:
                    try:
                        # TODO consider best way to override default func form
                        # On per bond, or per type basis
                        func_form = functional_forms[atomlist[-1]]
                    except AttributeError:
                        func_form = self._functional_forms[len(atomlist)]

                    if {x for x in atomlist if atomlist.count(x) > 1}:
                        raise ValueError("Defined bond '{0}' contains duplicate atoms".format(atomlist))

                    mol_bonds.append(Bond(atoms=atomlist, func_form=func_form))
                    if len(atomlist) > 2:
                        angles_defined = True

                if not angles_defined:
                    angles, dihedrals = self._create_angles(mol_bonds)

                    if options.generate_angles:
                        for atomlist in angles:
                            mol_bonds.append(Bond(atoms=atomlist, func_form=self._functional_forms[3]))

                    if options.generate_dihedrals:
                        for atomlist in dihedrals:
                            mol_bonds.append(Bond(atoms=atomlist, func_form=self._functional_forms[4]))

    @staticmethod
    def _create_angles(mol_bonds):
        """
        Create angles and dihedrals from bonded topology.

        :param mol_bonds: List of bonds within a molecule to generate angles for
        :return: List of angle and dihedral bond name tuples
        """
        bonds = [bond.atoms for bond in mol_bonds if len(bond.atoms) == 2]
        angles = extend_graph_chain(bonds, bonds)
        dihedrals = extend_graph_chain(angles, bonds)
        return angles, dihedrals

    def get_bond_lengths(self, mol, with_constr=False):
        """
        Return list of all bond lengths in molecule.  May include constraints.

        :param mol: Molecule name to return bonds for
        :param with_constr: Include constraints?
        :return: List of bonds
        """
        if with_constr:
            return [bond for bond in self._molecules[mol] if len(bond.atoms) == 2]
        else:
            return [bond for bond in self._molecules[mol] if len(bond.atoms) == 2 and bond.fconst < self._fconst_constr_threshold]

    def get_bond_length_constraints(self, mol):
        """
        Return list of all bond length constraints in molecule.

        :param mol: Molecule name to return bonds for
        :return: List of bonds
        """
        return [bond for bond in self._molecules[mol] if len(bond.atoms) == 2 and bond.fconst >= self._fconst_constr_threshold]

    def get_bond_angles(self, mol, exclude_triangle=True):
        """
        Return list of all bond angles in molecule.

        :param mol: Molecule name to return bonds for
        :param exclude_triangle: Exclude angles that are part of a triangle?
        :return: List of bonds
        """
        angles = [bond for bond in self._molecules[mol] if len(bond.atoms) == 3]

        if exclude_triangle:
            edges = [tuple(bond.atoms) for bond in self.get_bond_lengths(mol, with_constr=True)]

            def is_triangle(atoms):
                triangle_edges = 0
                for j in range(3):
                    if (atoms[j - 1], atoms[j]) in edges or (atoms[j], atoms[j - 1]) in edges:
                        triangle_edges += 1
                return triangle_edges >= 3

            angles = [angle for angle in angles if not is_triangle(angle.atoms)]

        return angles

    def get_bond_dihedrals(self, mol):
        """
        Return list of all bond dihedrals in molecule.

        :param mol: Molecule name to return bonds for
        :return: List of bonds
        """
        return [bond for bond in self._molecules[mol] if len(bond.atoms) == 4]

    def _populate_atom_numbers(self, mapping):
        """
        Add atom numbers to all bonds.

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
                    bond.atom_numbers = [index.index(atom.lstrip("+-")) for atom in bond.atoms]
                except ValueError as e:
                    missing = [atom for atom in bond.atoms if atom.lstrip("+-") not in index]
                    e.args = ("Bead(s) {0} do(es) not exist in residue {1}".format(missing, mol),)
                    raise

    def write_itp(self, filename, mapping):
        """
        Output a GROMACS .itp file containing atoms/beads and bonded terms.

        :param filename: Name of output file
        :param mapping: AA->CG Mapping from which to collect bead properties
        """
        self._populate_atom_numbers(mapping)
        backup_file(filename)

        def write_bond_angle_dih(bonds, section_header, itp, print_fconst=True, multiplicity=None, rad2deg=False):
            if bonds:
                print("\n[ {0:s} ]".format(section_header), file=itp)
            for bond in bonds:
                line = " ".join(["{0:4d}".format(atnum + 1) for atnum in bond.atom_numbers])
                eqm = math.degrees(bond.eqm) if rad2deg else bond.eqm
                line += " {0:4d} {1:12.5f}".format(bond.gromacs_type_id, eqm)
                if print_fconst:
                    line += " {0:12.5f}".format(bond.fconst)
                if multiplicity is not None:
                    line += " {0:4d}".format(multiplicity)
                print(line, file=itp)

        with open(filename, "w") as itp:
            header = ("; \n"
                      "; Topology prepared automatically using PyCGTOOL \n"
                      "; James Graham <J.A.Graham@soton.ac.uk> 2016 \n"
                      "; University of Southampton \n"
                      "; https://github.com/jag1g13/pycgtool \n"
                      ";")
            print(header, file=itp)

            # Print molecule
            not_calc = "  Parameters have not been calculated."
            for mol in self._molecules:
                if mol not in mapping:
                    logger.warning("Molecule '{0}' present in bonding file, but not in mapping.".format(mol) + not_calc)
                    continue
                if not all((bond.fconst is not None for bond in self._molecules[mol])):
                    logger.warning("Molecule '{0}' has no measured bond values.".format(mol) + not_calc)
                    continue

                print("\n[ moleculetype ]", file=itp)
                print("{0:4s} {1:4d}".format(mol, 1), file=itp)

                print("\n[ atoms ]", file=itp)
                for i, bead in enumerate(mapping[mol]):
                    #      atnum  type  resnum resname atname c-group  charge (mass)
                    print("{0:4d} {1:4s} {2:4d} {3:4s} {4:4s} {5:4d} {6:8.3f}".format(
                          i + 1, bead.type, 1, mol, bead.name, i + 1, bead.charge
                          ), file=itp)

                write_bond_angle_dih(self.get_bond_lengths(mol), "bonds", itp)
                write_bond_angle_dih(self.get_bond_angles(mol), "angles", itp, rad2deg=True)
                write_bond_angle_dih(self.get_bond_dihedrals(mol), "dihedrals", itp, multiplicity=1, rad2deg=True)
                write_bond_angle_dih(self.get_bond_length_constraints(mol), "constraints", itp, print_fconst=False)

    def apply(self, frame):
        """
        Calculate bond lengths/angles for a given Frame and store into Bonds.

        :param frame: Frame from which to calculate values
        """
        def calc_length(atoms):
            vec = dist_with_pbc(atoms[0].coords, atoms[1].coords, frame.box)
            return vector_len(vec)

        def calc_angle(atoms):
            veca = dist_with_pbc(atoms[0].coords, atoms[1].coords, frame.box)
            vecb = dist_with_pbc(atoms[1].coords, atoms[2].coords, frame.box)
            return math.pi - vector_angle(veca, vecb)

        def calc_dihedral(atoms):
            veca = dist_with_pbc(atoms[0].coords, atoms[1].coords, frame.box)
            vecb = dist_with_pbc(atoms[1].coords, atoms[2].coords, frame.box)
            vecc = dist_with_pbc(atoms[2].coords, atoms[3].coords, frame.box)

            c1 = vector_cross(veca, vecb)
            c2 = vector_cross(vecb, vecc)

            return vector_angle_signed(c1, c2, vecb)

        calc = {2: calc_length,
                3: calc_angle,
                4: calc_dihedral}

        for prev_res, res, next_res in sliding(frame):
            try:
                mol_meas = self._molecules[res.name]
            except KeyError:
                # Bonds have not been specified for molecule - probably water - ignore this residue
                continue

            adj_res = {"-": prev_res,
                       "+": next_res}

            for bond in mol_meas:
                try:
                    # TODO tidy this
                    atoms = [adj_res.get(name[0], res)[name.lstrip("-+")] for name in bond.atoms]
                    val = calc[len(atoms)](atoms)
                    bond.values.append(val)
                except (NotImplementedError, TypeError):
                    # TypeError is raised when residues on end of chain calc bond to next
                    pass
                except ZeroDivisionError as e:
                    e.args = ("Zero division in calculation of <{0}>".format(" ".join(bond.atoms)),)
                    raise e

    def boltzmann_invert(self, progress=False):
        """
        Perform Boltzmann Inversion of all bonds to calculate equilibrium value and force constant.

        :param progress: Display a progress bar using tqdm if available
        """
        bond_iter = itertools.chain(*self._molecules.values())
        bond_iter_wrap = bond_iter
        if progress:
            total = sum(map(len, self._molecules.values()))
            bond_iter_wrap = tqdm(bond_iter, total=total, ncols=80)

        for bond in bond_iter_wrap:
            try:
                bond.boltzmann_invert(temp=self._temperature)
            except ValueError:
                pass

    # TODO add test
    def dump_values(self, target_number=10000):
        """
        Output measured bond values to files for length, angles and dihedrals.

        :param target_number: Approx number of sample measurements to output.  If None, all samples will be output
        """

        def write_bonds_to_file(bonds, filename, rad2deg=False):
            with open(filename, "w") as f:
                for row in transpose_and_sample((bond.values for bond in bonds), n=target_number):
                    if rad2deg:
                        row = [math.degrees(val) for val in row]
                    print((len(row) * "{:12.5f}").format(*row), file=f)

        for mol in self._molecules:
            if mol == "SOL":
                continue
            bonds = self.get_bond_lengths(mol, with_constr=True)
            if bonds:
                write_bonds_to_file(bonds, "{0}_length.dat".format(mol))

            bonds = self.get_bond_angles(mol)
            if bonds:
                write_bonds_to_file(bonds, "{0}_angle.dat".format(mol), rad2deg=True)

            bonds = self.get_bond_dihedrals(mol)
            if bonds:
                write_bonds_to_file(bonds, "{0}_dihedral.dat".format(mol), rad2deg=True)

    def __len__(self):
        return len(self._molecules)

    def __contains__(self, item):
        return item in self._molecules

    def __getitem__(self, item):
        return self._molecules[item]

    def __iter__(self):
        return iter(self._molecules)
