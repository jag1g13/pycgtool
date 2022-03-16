"""
This module contains a single class ForceField used to output a GROMACS .ff forcefield.
"""

import os
import pathlib
import shutil
import typing

from .parsers import CFG
from .util import any_starts_with, backup_file, file_write_lines

PathLike = typing.Union[pathlib.Path, str]


def copy_files(
    src_dir: pathlib.Path, dest_dir: pathlib.Path, files: typing.Iterable[str]
) -> None:
    """Copy files from one directory to another."""
    for f in files:
        src_path = src_dir.joinpath(f)
        dest_path = dest_dir.joinpath(f)
        shutil.copyfile(str(src_path), str(dest_path))


class ForceField:
    """
    Class used to output a GROMACS .ff forcefield
    """

    def __init__(self, name: str, dir_path: PathLike = pathlib.Path(".")):
        """
        Open a named forcefield directory.  If it does not exist it is created.

        :param str name: Forcefield name to open/create
        """
        self.directory = pathlib.Path(dir_path).joinpath(f"ff{name}.ff")
        backup_file(self.directory)
        self.directory.mkdir(parents=True, exist_ok=True)

        with open(self.directory.joinpath("forcefield.itp"), "w") as itp:
            print(f"#define _FF_PYCGTOOL_{name}", file=itp)
            print('#include "martini_v2.2.itp"', file=itp)

        data_dir = pathlib.Path(__file__).parent.joinpath("data")

        # Copy MARTINI files
        copy_files(
            data_dir,
            self.directory,
            [
                "martini_v2.2.itp",
                "watermodels.dat",
                "w.itp",
            ],
        )

        # Create atomtypes.atp required for correct masses with pdb2gmx
        atomtypes_atp = os.path.join(self.directory, "atomtypes.atp")
        with CFG(data_dir.joinpath("martini_v2.2.itp")) as itp, open(
            atomtypes_atp, "w"
        ) as atomtypes:
            for toks in itp["atomtypes"]:
                print(" ".join(toks), file=atomtypes)

        with open(self.directory.joinpath("forcefield.doc"), "w") as doc:
            print(f"PyCGTOOL produced MARTINI force field - {name}", file=doc)

    def write(self, filename: str, mapping, bonds):
        """
        Write RTP and R2B files for this forcefield.

        :param str filename: Filename prefix to use for both files
        :param Mapping mapping: CG Mapping object
        :param Iterable[Bond] bonds: CG Bonds object
        """
        lines, nterms, cterms = ForceField.write_rtp(mapping, bonds)
        file_write_lines(self.directory.joinpath(f"{filename}.rtp"), lines)

        lines = ForceField.write_r2b(nterms, cterms)
        file_write_lines(self.directory.joinpath(f"{filename}.r2b"), lines)

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
    def needs_terminal_entries(mol_name, bondset):
        bonds = bondset.get_bonds(mol_name, natoms=-1)
        return any_starts_with(bonds, "-"), any_starts_with(bonds, "+")

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
        """

        def write_residue(mol_name, mol_mapping, strip=None, prepend=""):
            ret_lines = ["[ {0} ]".format(prepend + mol_name), "  [ atoms ]"]

            for bead in mol_mapping:
                #                      name   type  charge  chg-group
                ret_lines.append(
                    "    {:>4s} {:>4s} {:3.6f} {:4d}".format(
                        bead.name, bead.type, bead.charge, 0
                    )
                )

            for natoms, (section, multiplicity) in enumerate(
                (("bonds", None), ("angles", None), ("dihedrals", 1)), start=2
            ):
                if strip is None:
                    bond_list = bonds.get_bonds(mol_name, natoms)
                else:
                    bond_list = bonds.get_bonds(
                        mol_name,
                        natoms,
                        select=lambda bond: not any_starts_with(bond, strip),
                    )

                ret_lines.extend(
                    ForceField.bond_section(bond_list, section, multiplicity)
                )

            return ret_lines

        n_terms = set()
        c_terms = set()

        rtp_lines = ["[ bondedtypes ]", ("{:4d}" * 8).format(1, 1, 1, 1, 1, 1, 0, 0)]

        for mol_name, mol_mapping in mapping.items():
            try:
                rtp_lines.extend(write_residue(mol_name, mol_mapping))
            except KeyError:
                continue

            needs_terminal_entry = ForceField.needs_terminal_entries(mol_name, bonds)

            if needs_terminal_entry[0]:
                rtp_lines.extend(
                    write_residue(mol_name, mol_mapping, strip="-", prepend="N")
                )
                n_terms.add(mol_name)

            if needs_terminal_entry[1]:
                rtp_lines.extend(
                    write_residue(mol_name, mol_mapping, strip="+", prepend="C")
                )
                c_terms.add(mol_name)
                if needs_terminal_entry[0]:
                    rtp_lines.extend(
                        write_residue(
                            mol_name, mol_mapping, strip=("-", "+"), prepend="2"
                        )
                    )

        return rtp_lines, n_terms, c_terms

    @staticmethod
    def write_r2b(n_terms, c_terms):
        """
        Return lines of a GROMACS R2B file.

        This file defines names used for chain terminal records for PDB2GMX.

        :param Iterable[str] n_terms: Set of molecule names requiring N terminal records
        :param Iterable[str] c_terms: Set of molecule names requiring C terminal records
        :return List[str]: Lines of R2B file
        """
        ret_lines = [
            "; rtp residue to rtp building block table",
            ";     main  N-ter C-ter 2-ter",
        ]

        for resname in sorted(n_terms | c_terms):
            nter_str = ("N" + resname) if resname in n_terms else "-"
            cter_str = ("C" + resname) if resname in c_terms else "-"
            both_ter_str = ("2" + resname) if resname in (n_terms & c_terms) else "-"
            ret_lines.append(
                "{0:5s} {0:5s} {1:5s} {2:5s} {3:5s}".format(
                    resname, nter_str, cter_str, both_ter_str
                )
            )

        return ret_lines
