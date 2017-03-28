"""
This module contains a single class ForceField used to output a GROMACS .ff forcefield.
"""

import os
import shutil
import functools

from .util import dir_up, any_starts_with, file_write_lines
from .parsers import ITP


class ForceField:
    """
    Class used to output a GROMACS .ff forcefield
    """
    def __init__(self, name):
        """
        Open a named forcefield directory.  If it does not exist it is created.

        :param str name: Forcefield name to open/create
        """
        self.dirname = "ff{0}.ff".format(name)
        os.makedirs(self.dirname, exist_ok=True)

        with open(os.path.join(self.dirname, "forcefield.itp"), "w") as itp:
            print("#define _FF_PYCGTOOL_{0}".format(name), file=itp)
            print('#include "martini_v2.2.itp"', file=itp)

        data_dir = os.path.join(dir_up(os.path.realpath(__file__), 2), "data")
        # Copy main MARTINI itp
        martini_itp = os.path.join(data_dir, "martini_v2.2.itp")
        shutil.copyfile(martini_itp, os.path.join(self.dirname, "martini_v2.2.itp"))
        # Copy water models
        shutil.copyfile(os.path.join(data_dir, "watermodels.dat"),
                        os.path.join(self.dirname, "watermodels.dat"))
        shutil.copyfile(os.path.join(data_dir, "w.itp"),
                        os.path.join(self.dirname, "w.itp"))

        # Create atomtypes.atp required for correct masses with pdb2gmx
        atomtypes_atp = os.path.join(self.dirname, "atomtypes.atp")
        with ITP(martini_itp) as itp, open(atomtypes_atp, 'w') as atomtypes:
            for toks in itp["atomtypes"]:
                print(" ".join(toks), file=atomtypes)

        with open(os.path.join(self.dirname, "forcefield.doc"), "w") as doc:
            print("PyCGTOOL produced MARTINI force field - {0}".format(name), file=doc)

    def write(self, filename, mapping, bonds):
        """
        Write RTP and R2B files for this forcefield.

        :param str filename: Filename prefix to use for both files
        :param Mapping mapping: CG Mapping object
        :param Iterable[Bond] bonds: CG Bonds object
        """
        lines, nterms, cterms, bothterms = ForceField.write_rtp(mapping, bonds)
        file_write_lines(os.path.join(self.dirname, filename + ".rtp"), lines)

        lines = ForceField.write_r2b(nterms, cterms, bothterms)
        file_write_lines(os.path.join(self.dirname, filename + ".r2b"), lines)

    @staticmethod
    def bond_section(bonds, section_header, multiplicity=None):
        """
        Populate an RTP bond/angle/dihedral section.

        :param iterable[Bond] bonds: Iterable of bonds to add to RTP
        :param str section_header: RTP section header i.e. "bonds"/"angles"/"dihedrals"
        :param int multiplicity: Multiplicity of dihedral, default is None
        :return List[str]: Lines to add to RTP file
        """
        ret_lines = []
        if bonds:
            ret_lines.append("  [ {0:s} ]".format(section_header))
        for bond in bonds:
            line = "    " + " ".join(["{0:>4s}".format(atom) for atom in bond.atoms])
            line += " {0:12.5f} {1:12.5f}".format(bond.eqm, bond.fconst)
            if multiplicity is not None:
                line += " {0:4d}".format(multiplicity)
            ret_lines.append(line)
        return ret_lines

    @staticmethod
    def write_rtp(mapping, bonds):
        """
        Return lines of a GROMACS RTP file.

        This file defines the residues present in the forcefield and allows pdb2gmx to be used.

        :param Mapping mapping: AA->CG mapping from which to collect molecules
        :param BondSet bonds: BondSet from which to collect bonds
        :return (list[str], set[str], set[str], set[str]):
                List of lines for RTP file,
                Set of residues requiring N terminal records,
                Set of residues requiring C terminal records,
                Set of residues requiring N and C terminal records
        """

        def write_residue(name, strip=None, prepend=""):
            ret_lines = [
                "[ {0} ]".format(prepend + name),
                "  [ atoms ]"
            ]

            for bead in mapping[name]:
                #                      name   type  charge  chg-group
                ret_lines.append("    {:>4s} {:>4s} {:3.6f} {:4d}".format(
                    bead.name, bead.type, bead.charge, 0
                ))

            needs_terminal_entry = [False, False]

            get_bond_functions = [("bonds", functools.partial(bonds.get_bond_lengths, with_constr=True)),
                                  ("angles", bonds.get_bond_angles),
                                  ("dihedrals", bonds.get_bond_dihedrals)]

            for get_bond in get_bond_functions:
                bond_tmp = get_bond[1](name)
                if strip is not None:
                    bond_tmp = [bond for bond in bond_tmp if not any_starts_with(bond, strip)]
                ret_lines.extend(ForceField.bond_section(bond_tmp, get_bond[0],
                                                         multiplicity=1 if get_bond[0] == "dihedrals" else None))
                needs_terminal_entry[0] |= any_starts_with(bond_tmp, "-")
                needs_terminal_entry[1] |= any_starts_with(bond_tmp, "+")

            return ret_lines, needs_terminal_entry

        n_terms = set()
        c_terms = set()
        both_terms = set()

        rtp_lines = [
            "[ bondedtypes ]",
            ("{:4d}" * 8).format(1, 1, 1, 1, 1, 1, 0, 0)
        ]

        for mol in mapping:
            # Skip molecules not listed in bonds
            if mol not in bonds:
                continue

            needs_terminal_entry = write_residue(mol)
            if needs_terminal_entry[0]:
                rtp_lines.extend(write_residue(mol, strip="-", prepend="N")[0])
                n_terms.add(mol)
            if needs_terminal_entry[1]:
                rtp_lines.extend(write_residue(mol, strip="+", prepend="C")[0])
                c_terms.add(mol)
            if all(needs_terminal_entry):
                rtp_lines.extend(write_residue(mol, strip=("-", "+"), prepend="2")[0])
                both_terms.add(mol)

        return rtp_lines, n_terms, c_terms, both_terms

    @staticmethod
    def write_r2b(n_terms, c_terms, both_terms):
        """
        Return lines of a GROMACS R2B file.

        This file defines names used for chain terminal records for PDB2GMX.

        :param Iterable[str] n_terms: Set of molecule names requiring N terminal records
        :param Iterable[str] c_terms: Set of molecule names requiring C terminal records
        :param Iterable[str] both_terms: Set of molecule names requiring both terminal records
        :return List[str]: Lines of R2B file
        """
        ret_lines = [
            "; rtp residue to rtp building block table",
            ";     main  N-ter C-ter 2-ter"
        ]

        for resname in sorted(set.union(n_terms, c_terms, both_terms)):
            nter_str = ("N" + resname) if resname in n_terms else "-"
            cter_str = ("C" + resname) if resname in c_terms else "-"
            both_ter_str = ("2" + resname) if resname in both_terms else "-"
            ret_lines.append("{0:5s} {0:5s} {1:5s} {2:5s} {3:5s}".format(resname, nter_str, cter_str, both_ter_str))

        return ret_lines
