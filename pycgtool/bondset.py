"""
Module containing classes to calculate bonded properties from a Frame.

BondSet contains a dictionary of lists of Bonds.  Each list corresponds to a single molecule.
"""

import itertools
import math
import logging
import numpy as np
from copy import deepcopy


try:
    from tqdm import tqdm
except ImportError:
    from .util import tqdm_dummy as tqdm

from .mapping import VirtualMap
from .functionalforms import FunctionalForms
from .parsers.cfg import CFG
from .util import (
    circular_mean,
    circular_variance,
    dist_with_pbc,
    extend_graph_chain,
    file_write_lines,
    sliding,
    transpose_and_sample,
    vector_angle,
    vector_angle_signed,
    vector_cross,
    vector_len,
    merge_list_of_lists
)

logger = logging.getLogger(__name__)


class Molecule:
    """
    Holds data for a molecule comprised of multiple residues
    """
    __slots__ = ["resnames", "bonds"]

    def __init__(self, resnames=None, bonds=None):
        """
        :param resnames: list of residue names
        :param bonds:  list of lists of Bond objects in the same order as 'resnames'
        """
        self.resnames = resnames
        self.bonds = bonds
        if bonds is None:
            self.bonds = []


    def add_bond(self, bond):
        """
        Add bond to object

        :param bond: instance of class with 'Bond' as their base class
        """
        self.bonds.append(bond)

    @property
    def inter(self):
        """collects inter residue bonds"""
        inter = []
        for bond in self.bonds:
            if isinstance(bond, GlobalBond):
                inter.append(bond)
        return inter

    @property
    def is_multiresidue(self):
        """checks if Molecule has multiple residues"""
        return True if len(self.inter) > 0 else False

    def __len__(self):
        return len(self.bonds)

    def __iter__(self):
        return iter(self.bonds)

    def __getitem__(self, item):
        return self.bonds[item]

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
            # if dihedral get value to shift cosine function by, NOT equilibrium value
            if len(self.atoms) == 4:
                two_pi = 2 * np.pi
                self.eqm -= np.pi
                self.eqm = (((self.eqm + np.pi) % two_pi) - np.pi)


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

class GlobalBond(Bond):
    __slots__ = ["atoms", "atom_numbers", "values", "eqm", "fconst", "gromacs_type_id", "_func_form",
                 "resids", "resnames", "labels"]
    def __init__(self, atoms, atom_numbers=None, func_form=None):
        """
        Class for bonds between not necessarily adjacent residues
        :param List[str] atoms: List of atom names defining the bond
        :param List[int] atom_numbers: List of atom numbers defining the bond
        :param func_form: Functional form to use for Boltzmann Inversion
        """
        self.labels = atoms
        atom_names = []
        resids = []
        resnames = []
        for atom in atoms:
            try:
                name, resid, resname = atom.split("_")
                resid = int(resid)
                atom_names.append(name)
                resids.append(resid)
                resnames.append(resname)

            except ValueError:
                print("incorrect syntax for entry '{}' in [global] section".format(atom))
                raise SyntaxError
        Bond.__init__(self, atom_names, atom_numbers=atom_numbers, func_form=func_form)
        self.resids = resids
        self.resnames = resnames

    def get_atoms(self, frame):
        """
        get atoms involved in bond
        :param frame: Frame object
        :return: list of Atom objects
        """
        atoms = []
        for resid, resname, name in zip(self.resids, self.resnames, self.atoms):
            for residue in frame:
                if residue.num == resid:
                    if residue.name == resname:
                        try:
                            atom = residue[name]
                            atoms.append(atom)
                        except KeyError:
                            pass
        return atoms

    def get_residue_ids(self, frame):
        """
        Get internal resids of residues involved with bond
        :param frame:
        :return: list of resids
        """
        residue_ids = []
        for resid, resname, name in zip(self.resids, self.resnames, self.atoms):
            for ind, residue in enumerate(frame):
                if residue.num == resid:
                    if residue.name == resname:
                        residue_ids.append(ind)
        return residue_ids

    def populate_ids(self, mol_beads):
        """
        populate internal indices of bond
        :param mol_beads: beads in the molecule
        """
        ids = []
        for resid, resname, name in zip(self.resids, self.resnames, self.atoms):
            for ind, beads in enumerate(mol_beads, start=1):
                index = [bead.name for bead in beads]
                if ind == resid:
                    try:
                       ids.extend([ beads[index.index(atom)].num for atom in index if atom == name])
                    except ValueError as e:
                        missing = [atom for atom in self.atoms if atom.lstrip("+-") not in index]
                        e.args = ("Bead(s) {0} do(es) not exist in residue {1}".format(missing, resname),)
                        raise
                    except KeyError:
                        pass
        self.atom_numbers = ids

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
        functional_forms_map = FunctionalForms()
        if self._default_fc:
            default_forms = ["MartiniDefaultLength", "MartiniDefaultAngle", "MartiniDefaultDihedral"]
        else:
            default_forms = ["Harmonic", "CosHarmonic", "Harmonic"]
        self._functional_forms = [None, None]
        self._functional_forms.extend(map(lambda x: functional_forms_map[x], default_forms))

        try:
            self._functional_forms[2] = functional_forms_map[options.length_form]
        except AttributeError:
            pass

        try:
            self._functional_forms[3] = functional_forms_map[options.angle_form]
        except AttributeError:
            pass

        try:
            self._functional_forms[4] = functional_forms_map[options.dihedral_form]
        except AttributeError:
            pass

        molecule_delimiters = ('<', '>')
        is_global = False
        multi_count = 0
        with CFG(filename) as cfg:
            for mol_name, mol_section in cfg.items():
                if "<" in mol_name:
                    is_global = True
            for mol_name, mol_section in cfg.items():
                now_global = False
                if "<" in mol_name:
                    multi_count += 1
                    if multi_count > 1:
                        print("More than one multi residue molecule detected - this is not yet supported by pycgtool")
                        raise NotImplementedError
                    now_global = True
                    mol_name = mol_name.strip(" ".join(molecule_delimiters))
                    self._molecules[mol_name] = Molecule()

                    mol_bonds = self._molecules[mol_name]

                else:
                    self._molecules[mol_name] = Molecule(resnames=[mol_name])
                    mol_bonds = self._molecules[mol_name]

                angles_defined = False
                for atomlist in mol_section:
                    is_angle_or_dihedral = len(atomlist) > 2
                    mean_function = circular_mean if is_angle_or_dihedral else np.nanmean
                    variance_function = circular_variance if is_angle_or_dihedral else np.nanvar

                    # Construct instance of Boltzmann Inversion function and
                    # inject dependencies for mean and variance functions
                    try:
                        # TODO consider best way to override default func form
                        # On per bond, or per type basis
                        func_form = functional_forms_map[atomlist[-1]](
                            mean_function,
                            variance_function
                        )
                    except AttributeError:
                        func_form = self._functional_forms[len(atomlist)](
                            mean_function,
                            variance_function
                        )

                    if {x for x in atomlist if atomlist.count(x) > 1}:
                        raise ValueError("Defined bond '{0}' contains duplicate atoms".format(atomlist))
                    if now_global:
                        mol_bonds.add_bond(GlobalBond(atoms=atomlist, func_form=func_form))
                    else:
                        mol_bonds.add_bond(Bond(atoms=atomlist, func_form=func_form))
                    if len(atomlist) > 2:
                        angles_defined = True

                if not angles_defined:
                    if is_global:
                        if options.generate_angles or options.generate_dihedrals:
                            logger.warning("Automated generation of angles or dihedrals between "
                                           "residues not implemented! Please specifiy angles and dihedrals in [< {0} >]"
                                           "section of *.bnd file".format(mol_name))
                    else:

                        angles, dihedrals = self._create_angles(mol_bonds)

                        if options.generate_angles:
                            for atomlist in angles:
                                mol_bonds.add_bond(
                                    Bond(
                                        atoms=atomlist,
                                        func_form=self._functional_forms[3](circular_mean, circular_variance)
                                    )
                                )

                        if options.generate_dihedrals:
                            for atomlist in dihedrals:
                                mol_bonds.add_bond(
                                    Bond(
                                        atoms=atomlist,
                                        func_form=self._functional_forms[4](circular_mean, circular_variance)
                                    )
                                )

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

    def get_bonds(self, mol, natoms, select=lambda x: True):
        """
        Return list of bonds from molecule containing natoms atoms.

        :param str mol: Molecule name
        :param int natoms: Number of atoms in bond, i.e. 2 for bond length, 3 for angle, 4 for dihedral
        :param function select: Optional lambda, return only bonds for which this is True
        :return List[Bond]: List of bonds
        """
        bonds = self._molecules[mol]

        if natoms == -1:
            return [bond for bond in bonds if select(bond)]
        return [bond for bond in bonds if len(bond.atoms) == natoms and select(bond)]

    def get_bond_lengths(self, mol, with_constr=False):
        """
        Return list of all bond lengths in molecule.  May include constraints.

        :param mol: Molecule name to return bonds for
        :param with_constr: Include constraints?
        :return: List of bonds
        """
        bonds = self._molecules[mol]

        if with_constr:
            return [bond for bond in bonds if len(bond.atoms) == 2]
        else:
            return [bond for bond in bonds if len(bond.atoms) == 2 and bond.fconst < self._fconst_constr_threshold]

    def get_bond_length_constraints(self, mol):
        """
        Return list of all bond length constraints in molecule.

        :param mol: Molecule name to return bonds for
        :return: List of bonds
        """
        bonds = self._molecules[mol]
        return [bond for bond in bonds if len(bond.atoms) == 2 and bond.fconst >= self._fconst_constr_threshold]

    def get_bond_angles(self, mol, exclude_triangle=True):
        """
        Return list of all bond angles in molecule.

        :param mol: Molecule name to return bonds for
        :param exclude_triangle: Exclude angles that are part of a triangle?
        :return: List of bonds
        """
        bonds = self._molecules[mol]
        angles = [bond for bond in bonds if len(bond.atoms) == 3]

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
        bonds = self._molecules[mol]
        return [bond for bond in bonds if len(bond.atoms) == 4]

    def get_virtual_beads(self, mapping):
        """
        Return list of all virtual beads in molecule
        :param mapping:
        :return: list of virtual beads
        """
        return [bead for bead in mapping if isinstance(bead, VirtualMap)]

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

    def connect_residues(self, frame, mapping):
        """
        connects residues together to form new molecules
        :param frame: Frame from which to determine inter residue connections
        """
        # TODO this section is a bit verbose - simplify
        not_multi = [mol for mol in self._molecules if not self._molecules[mol].is_multiresidue]
        for mol in self._molecules:
            if self._molecules[mol].is_multiresidue:
                mol_bonds = self._molecules[mol].inter
                fragments_resid = []
                for bond in mol_bonds:
                    fragments_resid.append(bond.get_residue_ids(frame))

                self._populate_atom_numbers(mapping)

                molecule_internal_resids = merge_list_of_lists(fragments_resid)
                if len(molecule_internal_resids) > 1:
                    print("All fragments of Molecule '{}' are not connected - add missing bonds to .bnd".format(mol))
                    raise SyntaxError

                # populate residue bond ids
                molecule_internal_resids = molecule_internal_resids[0]
                resnames = [frame[resid].name for resid in molecule_internal_resids]
                all_bonds = []
                all_beads = []
                start = 0
                for resid, resname in enumerate(resnames):
                    if resname in mapping:
                        beads = list(map(deepcopy, mapping[resname]))
                        if len(not_multi) > 0:
                            bonds = list(map(deepcopy, self._molecules[resname]))
                            all_bonds.extend(bonds)

                        for i, bead in enumerate(beads):
                            bead.num = start + i

                        if len(not_multi) > 0:
                            index = [bead.name for bead in beads]
                            for bond in bonds:
                                try:
                                    bond.atom_numbers = [index.index(atom.lstrip("+-")) for atom in bond.atoms]
                                except ValueError as e:
                                    missing = [atom for atom in bond.atoms if atom.lstrip("+-") not in index]
                                    e.args = ("Bead(s) {0} do(es) not exist in residue {1}".format(missing, resname),)
                                    raise
                                bond.atom_numbers = [start + ind for ind in bond.atom_numbers]

                        all_beads.append(beads)
                        start = beads[-1].num + 1

                # populate inter-residue bonds ids
                for bond in mol_bonds:
                    bond.populate_ids(all_beads)

                all_bonds.extend(mol_bonds)
                molecule = Molecule(resnames, all_bonds)
                self._molecules[mol] = molecule
        return

    def write_itp(self, filename, mapping):
        """
        Output a GROMACS .itp file containing atoms/beads and bonded terms.

        :param filename: Name of output file
        :param mapping: AA->CG Mapping from which to collect bead properties
        """
        self._populate_atom_numbers(mapping)
        file_write_lines(filename, self.itp_text(mapping))

    def itp_text(self, mapping):
        atom_template = {"nomass": "{0:4d} {1:4s} {2:4d} {3:4s} {4:4s} {5:4d} {6:8.3f}",
                         "mass": "{0:4d} {1:4s} {2:4d} {3:4s} {4:4s} {5:4d} {6:8.3f} {7:8.3f}"}

        def write_bond_angle_dih(bonds, section_header, print_fconst=True, multiplicity=None, rad2deg=False):
            ret_lines = []
            if bonds:
                ret_lines.append("\n[ {0:s} ]".format(section_header))
            for bond in bonds:
                line = " ".join(["{0:4d}".format(atnum + 1) for atnum in bond.atom_numbers])
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
        ]

        # Print molecule
        not_calc = "  Parameters have not been calculated."
        for mol in self._molecules:
            molecule = self._molecules[mol]

            if not all((bond.fconst is not None for bond in self._molecules[mol])):
                logger.warning("Molecule '{0}' has no measured bond values.".format(mol) + not_calc)
                continue

            ret_lines.append("\n[ moleculetype ]")
            ret_lines.append("{0:4s} {1:4d}".format(mol, 1))

            ret_lines.append("\n[ atoms ]")

            # print residues
            start = 1
            for resid, res in enumerate(molecule.resnames, start=1):
                beads = mapping[res]


                if res not in mapping:
                    logger.warning("Residue '{0}' present in bonding file, but not in mapping.".format(mol) + not_calc)
                    continue

                for i, bead in enumerate(beads, start=1):
                    #                 atnum   type resnum resname atname c-group charge (mass)
                    if isinstance(bead, VirtualMap):
                        ret_lines.append(atom_template["mass"].format(
                            start + bead.num, bead.type, resid, res, bead.name, start + bead.num, bead.charge, bead.mass
                        ))

                    else:
                        ret_lines.append(atom_template["nomass"].format(
                            start + bead.num, bead.type, resid, res, bead.name, start + bead.num, bead.charge
                        ))
                start = beads[-1].num + 2

            # print virtual sites
            virtual_beads = [self.get_virtual_beads(mapping[res]) for res in self._molecules[mol].resnames]
            virtual_beads = [bead for res_beads in virtual_beads for bead in res_beads]
            if len(virtual_beads) != 0:
                ret_lines.append("\n[ virtual_sitesn ]")
                excl_lines = ["\n[ exclusions ]"]
                for res in molecule.resnames:
                    beads = mapping[res]
                    virtual_beads = self.get_virtual_beads(beads)
                    for vbead in virtual_beads:
                        cg_ids = [bead.num + 1 for bead in beads if bead.name in vbead.atoms]
                        cg_ids.sort()
                        cg_ids_string = " ".join(map(str, cg_ids))
                        ret_lines.append(
                            "{0:^6d} {1:^6d} {2}".format(vbead.num + 1, vbead.gromacs_type_id, cg_ids_string))
                        vsite_exclusions = "{} ".format(vbead.num + 1) + cg_ids_string
                        excl_lines.append(vsite_exclusions)
                ret_lines.extend(excl_lines)

            ret_lines.extend(write_bond_angle_dih(self.get_bond_lengths(mol), "bonds"))
            ret_lines.extend(write_bond_angle_dih(self.get_bond_angles(mol), "angles", rad2deg=True))
            ret_lines.extend(write_bond_angle_dih(self.get_bond_dihedrals(mol), "dihedrals", multiplicity=1,
                                                  rad2deg=True))
            ret_lines.extend(write_bond_angle_dih(self.get_bond_length_constraints(mol), "constraints",
                                                  print_fconst=False))

        return ret_lines

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

        # TODO tidy this section
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
                    atoms = [adj_res.get(name[0], res)[name.lstrip("-+")] for name in bond.atoms]
                    num_atoms = len(atoms)
                    if num_atoms == 0:
                        raise Exception("Error: {0} residue in *.map and *.gro file, but atoms {1} not in gro".format(
                            res.name, " ".join(bond.atoms)))
                    val = calc[num_atoms](atoms)
                    bond.values.append(val)
                except (NotImplementedError, TypeError):
                    # TypeError is raised when residues on end of chain calc bond to next
                    pass
                except ZeroDivisionError as e:
                    e.args = ("Zero division in calculation of <{0}>".format(" ".join(bond.atoms)),)
                    raise e

        #calculate values for global bonds
        for mol in self._molecules:
            if self._molecules[mol].is_multiresidue:
                bonds = self._molecules[mol].bonds
                for bond in bonds:
                    try:
                        atoms = bond.get_atoms(frame)
                        num_atoms = len(atoms)
                        if num_atoms == 0:
                            raise Exception(
                                "Error: {0} Molecule in *.map, but atoms {1} not in gro".format(
                                    mol, " ".join(bond.atoms)))
                        val = calc[num_atoms](atoms)
                        bond.values.append(val)
                    except ZeroDivisionError as e:
                        e.args = ("Zero division in calculation of <{0}>".format(" ".join(bond.labels)),)
                        raise e

    def boltzmann_invert(self, progress=False):
        """
        Perform Boltzmann Inversion of all bonds to calculate equilibrium value and force constant.

        :param progress: Display a progress bar using tqdm if available
        """
        bond_iter = itertools.chain(*self._molecules.values())
        for mol in self._molecules:
            if self._molecules[mol].is_multiresidue:
                bonds = self._molecules[mol].bonds
                bond_iter = itertools.chain(bond_iter, bonds)
        if progress:
            total = sum(map(len, self._molecules.values()))
            bond_iter = tqdm(bond_iter, total=total, ncols=80)

        for bond in bond_iter:
            try:
                bond.boltzmann_invert(temp=self._temperature)
            except ValueError:
                pass

    @staticmethod
    def _get_lines_for_bond_dump(bonds, target_number=None, rad2deg=False):
        """
        Return a dump of bond measurements as a list of lines of text.

        :param bonds: Iterable of bonds
        :param int target_number: Sample size of bonds to dump - None will dump all bonds
        :param bool rad2deg: Should bond measurements be converted from radians to degrees?
        :return List[str]: Lines of bond measurements dump
        """
        ret_lines = []
        for row in transpose_and_sample((bond.values for bond in bonds), n=target_number):
            if rad2deg:
                row = [math.degrees(val) for val in row]
            ret_lines.append((len(row) * "{:12.5f}").format(*row))
        return ret_lines

    def dump_values(self, target_number=10000):
        """
        Output measured bond values to files for length, angles and dihedrals.

        :param int target_number: Approx number of sample measurements to output.  If None, all samples will be output
        """

        for mol in self._molecules:
            if mol == "SOL":
                continue
            if not self._molecules[mol].is_multiresidue:

                bonds = self.get_bond_lengths(mol, with_constr=True)
                if bonds:
                    lines = BondSet._get_lines_for_bond_dump(bonds, target_number)
                    file_write_lines("{0}_length.dat".format(mol), lines)

                bonds = self.get_bond_angles(mol)
                if bonds:
                    lines = BondSet._get_lines_for_bond_dump(bonds, target_number, rad2deg=True)
                    file_write_lines("{0}_angle.dat".format(mol), lines)

                bonds = self.get_bond_dihedrals(mol)
                if bonds:
                    lines = BondSet._get_lines_for_bond_dump(bonds, target_number, rad2deg=True)
                    file_write_lines("{0}_dihedral.dat".format(mol), lines)
            else:
                bonds = self.get_bond_lengths(mol, with_constr=True)
                global_bonds = [bond for bond in bonds if isinstance(bond, GlobalBond)]
                if global_bonds:
                    lines = BondSet._get_lines_for_bond_dump(global_bonds, target_number)
                    file_write_lines("{0}_length.dat".format(mol), lines)

                bonds = self.get_bond_angles(mol)
                global_bonds = [bond for bond in bonds if isinstance(bond, GlobalBond)]
                if global_bonds:
                    lines = BondSet._get_lines_for_bond_dump(bonds, target_number, rad2deg=True)
                    file_write_lines("{0}_angle.dat".format(mol), lines)

                bonds = self.get_bond_dihedrals(mol)
                global_bonds = [bond for bond in bonds if isinstance(bond, GlobalBond)]
                if global_bonds:
                    lines = BondSet._get_lines_for_bond_dump(bonds, target_number, rad2deg=True)
                    file_write_lines("{0}_dihedral.dat".format(mol), lines)



    def __len__(self):
        return len(self._molecules)

    def __contains__(self, item):
        return item in self._molecules

    def __getitem__(self, item):
        return self._molecules[item]

    def __iter__(self):
        return iter(self._molecules)
