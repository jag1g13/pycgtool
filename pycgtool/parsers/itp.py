"""Module containing classes used to parse ITP/TOP file format."""

import collections

from .cfg import CFG, NoSectionError, DuplicateSectionError


class ITP(CFG):
    """Class representing an .itp file

    Contains a dictionary for every molecule definition
    """
    def _read_file(self) -> None:
        mol_sections = ["atoms", "bonds", "angles", "dihedrals"]

        with open(self.filename) as f:
            curr_section = None
            curr_mol = None

            for line in f:
                line = self._read_line(line)
                if line is None:
                    continue

                if line.startswith("["):
                    curr_section = line.strip("[ ]")
                    if curr_section in mol_sections:
                        if curr_section not in self[curr_mol]:
                            self[curr_mol][curr_section] = []
                    continue

                toks = tuple(line.split())
                if curr_section == "moleculetype":
                    curr_mol = toks[0]
                    if curr_mol in self:
                        raise DuplicateSectionError(curr_mol, self.filename)
                    self[curr_mol] = collections.OrderedDict()

                elif curr_section in mol_sections:
                    try:
                        self[curr_mol][curr_section].append(toks)
                    except KeyError as e:
                        raise NoSectionError(self.filename) from e
