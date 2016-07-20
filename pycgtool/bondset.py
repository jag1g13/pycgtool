"""
Module containing classes to calculate bonded properties from a Frame.

BondSet contains a dictionary of lists of Bonds.  Each list corresponds to a single molecule.
"""

import itertools

import numpy as np
import math

try:
    from tqdm import tqdm
except ImportError:
    pass

from .util import stat_moments, sliding, dist_with_pbc, transpose_and_sample
from .util import extend_graph_chain, cross, backup_file
from .parsers.cfg import CFG

np.seterr(all="raise")


class Bond:
    """
    Class holding the properties of a single bonded term.

    Bond lengths, angles and dihedrals are all equivalent, distinguised by the number of atoms present.
    """
    __slots__ = ["atoms", "atom_numbers", "values", "eqm", "fconst"]

    def __init__(self, atoms=None, atom_numbers=None):
        """
        Create a single bond definition.

        :param atoms: List of atom names defining the bond
        :param atom_numbers: List of atom numbers defining the bond
        :return: Instance of Bond
        """
        self.atoms = atoms
        self.atom_numbers = atom_numbers
        self.values = []
        self.eqm = None
        self.fconst = None

    def __len__(self):
        return len(self.atoms)

    def boltzmann_invert(self, temp=310, angle_default_fc=True):
        """
        Perform Boltzmann Inversion using measured values of this bond to calculate equilibrium value and force constant.

        :param temp: Temperature at which the simulation was performed
        :param angle_default_fc: Use default value of 25 kJ mol-1 rad-2 for angle force constants
        """
        mean, var = stat_moments(self.values)

        rt = 8.314 * temp / 1000.
        rad2 = np.pi * np.pi / (180. * 180.)
        conv = {2: lambda: rt / var,
                3: lambda: 25. if angle_default_fc else rt * np.sin(np.radians(mean))**2 / (var * rad2),
                # 3: lambda: 25. if angle_default_fc else rt / (var * rad2),
                4: lambda: rt / (var * rad2)}

        self.eqm = mean
        try:
            self.fconst = conv[len(self.atoms)]()
        except FloatingPointError:
            self.fconst = 0

    def r_squared(self):
        raise NotImplementedError("Bond r-squared is not yet implemented")

    def __repr__(self):
        try:
            return "<Bond containing atoms {0} with r_0 {1:.3f} and force constant {2:.3e}>".format(
                ", ".join(self.atoms), self.eqm, self.fconst
            )
        except (AttributeError, TypeError):
            return "<Bond containing atoms {0}>".format(", ".join(self.atoms))


def angle(a, b):
    """
    Calculate the angle between two vectors.

    :param a: First vector
    :param b: Second vector
    :return: Angle in radians
    """
    dot = a[0]*b[0] + a[1]*b[1] + a[2]*b[2]
    mag = math.sqrt((a[0]*a[0] + a[1]*a[1] + a[2]*a[2]) * (b[0]*b[0] + b[1]*b[1] + b[2]*b[2]))
    try:
        return math.acos(dot / mag)
    except ValueError as e:
        val = dot / mag
        if abs(val) - 1 < 0.01:
            # If within 1% of acceptable value correct it, else reraise
            return math.acos(max(-1, min(1, dot / mag)))
        e.args = (e.args[0] + " in acos(" + str(dot / mag) + ")",)
        raise


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
            self._angle_default_fc = options.angle_default_fc
        except AttributeError:
            self._angle_default_fc = True

        with CFG(filename) as cfg:
            for mol in cfg:
                self._molecules[mol.name] = []
                molbnds = self._molecules[mol.name]

                angles_defined = False
                for atomlist in mol:
                    molbnds.append(Bond(atoms=atomlist))
                    if len(atomlist) > 2:
                        angles_defined = True

                if not angles_defined and options.generate_angles:
                    self._create_angles(mol.name, options.generate_dihedrals)

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

    def _create_angles(self, mol, gen_dihedrals=False):
        """
        Create angles (and dihedrals) from bonded topology.

        :param mol: Molecule to generate angles for
        :param gen_dihedrals: Also generate dihedrals?
        """
        bonds = [bond.atoms for bond in self._molecules[mol] if len(bond.atoms) == 2]

        angles = extend_graph_chain(bonds, bonds)
        for atomlist in angles:
            self._molecules[mol].append(Bond(atoms=atomlist))

        if gen_dihedrals:
            dihedrals = extend_graph_chain(angles, bonds)
            for atomlist in dihedrals:
                self._molecules[mol].append(Bond(atoms=atomlist))

    def _populate_atom_numbers(self, mapping):
        """
        Add atom numbers to all bonds.

        Uses previously defined atom names.

        :param mapping: AA->CG mapping used to collect atom/bead numbers
        """
        for mol in self._molecules:
            molmap = mapping[mol]
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
        backup_file(filename, verbose=True)

        with open(filename, "w") as itp:
            header = ("; \n"
                      "; Topology prepared automatically using PyCGTOOL \n"
                      "; James Graham <J.A.Graham@soton.ac.uk> 2016 \n"
                      "; University of Southampton \n"
                      "; https://github.com/jag1g13/pycgtool \n"
                      ";")
            print(header, file=itp)

            # Print molecule
            for mol in self._molecules:
                print("\n[ moleculetype ]", file=itp)
                print("{0:4s} {1:4d}".format(mol, 1), file=itp)

                print("\n[ atoms ]", file=itp)
                for i, bead in enumerate(mapping[mol]):
                    #      atnum  type  resnum resname atname c-group  charge (mass)
                    print("{0:4d} {1:4s} {2:4d} {3:4s} {4:4s} {5:4d} {6:8.3f}".format(
                          i + 1, bead.type, 1, mol, bead.name, i + 1, bead.charge
                          ), file=itp)

                bonds = self.get_bond_lengths(mol)
                if len(bonds):
                    print("\n[ bonds ]", file=itp)
                for bond in bonds:
                    print("{0:4d} {1:4d} {2:4d} {3:12.5f} {4:12.5f}".format(
                        bond.atom_numbers[0] + 1, bond.atom_numbers[1] + 1,
                        1, bond.eqm, bond.fconst
                    ), file=itp)

                bonds = self.get_bond_angles(mol)
                if len(bonds):
                    print("\n[ angles ]", file=itp)
                for bond in bonds:
                    print("{0:4d} {1:4d} {2:4d} {3:4d} {4:12.5f} {5:12.5f}".format(
                        bond.atom_numbers[0] + 1, bond.atom_numbers[1] + 1, bond.atom_numbers[2] + 1,
                        1, bond.eqm, bond.fconst
                    ), file=itp)

                bonds = self.get_bond_dihedrals(mol)
                if len(bonds):
                    print("\n[ dihedrals ]", file=itp)
                for bond in bonds:
                    print("{0:4d} {1:4d} {2:4d} {3:4d} {4:4d} {5:12.5f} {6:12.5f} {7:4d}".format(
                        bond.atom_numbers[0] + 1, bond.atom_numbers[1] + 1,
                        bond.atom_numbers[2] + 1, bond.atom_numbers[3] + 1,
                        1, bond.eqm, bond.fconst, 1
                    ), file=itp)

                bonds = self.get_bond_length_constraints(mol)
                if len(bonds):
                    print("\n[ constraints ]", file=itp)
                for bond in bonds:
                    print("{0:4d} {1:4d} {2:4d} {3:12.5f}".format(
                        bond.atom_numbers[0] + 1, bond.atom_numbers[1] + 1,
                        1, bond.eqm
                    ), file=itp)

    def apply(self, frame):
        """
        Calculate bond lengths/angles for a given Frame and store into Bonds.

        :param frame: Frame from which to calculate values
        """
        def calc_length(atoms):
            vec = dist_with_pbc(atoms[0].coords, atoms[1].coords, frame.box)
            return math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2])

        def calc_angle(atoms):
            veca = dist_with_pbc(atoms[0].coords, atoms[1].coords, frame.box)
            vecb = dist_with_pbc(atoms[1].coords, atoms[2].coords, frame.box)
            ret = np.degrees(np.pi - angle(veca, vecb))
            return ret

        def calc_dihedral(atoms):
            veca = dist_with_pbc(atoms[0].coords, atoms[1].coords, frame.box)
            vecb = dist_with_pbc(atoms[1].coords, atoms[2].coords, frame.box)
            vecc = dist_with_pbc(atoms[2].coords, atoms[3].coords, frame.box)

            c1 = cross(veca, vecb)
            c2 = cross(vecb, vecc)
            c3 = cross(c1, c2)

            ang = np.degrees(angle(c1, c2))
            direction = np.dot(vecb, c3)
            return ang if direction > 0 else -ang

        calc = {2: calc_length,
                3: calc_angle,
                4: calc_dihedral}

        for prev_res, res, next_res in sliding(frame):
            try:
                mol_meas = self._molecules[res.name]
            except KeyError:
                # Bonds have not been specified for molecule - probably water
                continue

            adj_res = {"-": prev_res,
                       "+": next_res}

            for bond in mol_meas:
                try:
                    # TODO tidy this
                    atoms = [adj_res.get(name[0], res)[name.lstrip("-+")] for name in bond.atoms]
                    val = calc[len(atoms)](atoms)
                    bond.values.append(val)
                except (NotImplementedError, TypeError, FloatingPointError):
                    # NotImplementedError is raised if form is not implemented
                    # TypeError is raised when residues on end of chain calc bond to next
                    pass

    def boltzmann_invert(self, progress=False):
        """
        Perform Boltzmann Inversion of all bonds to calculate equilibrium value and force constant.

        :param progress: Display a progress bar using tqdm if available
        """
        bond_iter = itertools.chain(*self._molecules.values())
        bond_iter_wrap = bond_iter
        if progress:
            try:
                total = sum(map(len, self._molecules.values()))
                bond_iter_wrap = tqdm(bond_iter, total=total)
            except NameError:
                pass

        for bond in bond_iter_wrap:
            bond.boltzmann_invert(temp=self._temperature,
                                  angle_default_fc=self._angle_default_fc)

    def dump_values(self, target_number=100000):
        """
        Output measured bond values to files for length, angles and dihedrals.

        :param target_number: Approx number of sample measurements to output.  If None, all samples will be output
        """

        for mol in self._molecules:
            if mol == "SOL":
                continue
            with open("{0}_length.dat".format(mol), "w") as f:
                bonds = self.get_bond_lengths(mol, with_constr=True)
                for row in transpose_and_sample((bond.values for bond in bonds), n=target_number):
                    print((len(row) * "{:12.5f}").format(*row), file=f)

            with open("{0}_angle.dat".format(mol), "w") as f:
                bonds = self.get_bond_angles(mol)
                for row in transpose_and_sample((bond.values for bond in bonds), n=target_number):
                    print((len(row) * "{:12.5f}").format(*row), file=f)

            with open("{0}_dihedral.dat".format(mol), "w") as f:
                bonds = self.get_bond_dihedrals(mol)
                for row in transpose_and_sample((bond.values for bond in bonds), n=target_number):
                    print((len(row) * "{:12.5f}").format(*row), file=f)

    def __len__(self):
        return len(self._molecules)

    def __contains__(self, item):
        return item in self._molecules

    def __getitem__(self, item):
        return self._molecules[item]

    def __iter__(self):
        return iter(self._molecules)
