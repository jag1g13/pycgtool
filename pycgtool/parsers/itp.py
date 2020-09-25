"""Module containing classes used to parse ITP/TOP file format."""

import collections
import logging

from .cfg import CFG, NoSectionError, DuplicateSectionError

logger = logging.getLogger(__name__)


class ITP(CFG):
    """Class representing an .itp file

    Contains a dictionary for every molecule definition
    """
    def _read_file(self) -> None:
        mol_sections = ["atoms", "bonds", "angles", "dihedrals"]

        with open(self.filename) as itp_file:
            curr_section = None
            curr_mol = None

            for line in itp_file:
                line = self._read_line(line)

                if line is None:
                    continue

                if line.startswith("["):
                    curr_section = line.strip("[ ]")

                    if curr_section in mol_sections:  # pragma: no cover
                        # This allows a section to appear twice within a single molecule definition
                        # It is not clear whether this is allowed by GROMACS but does no harm
                        self[curr_mol].setdefault(curr_section, [])

                    continue  # coverage.py cannot detect this branch due to the peephole optimizer

                toks = tuple(line.split())

                if curr_section == "moleculetype":
                    curr_mol = toks[0]

                    if curr_mol in self:
                        raise DuplicateSectionError(curr_mol, self.filename)

                    self[curr_mol] = collections.OrderedDict()

                elif curr_section in mol_sections:
                    self[curr_mol][curr_section].append(toks)

                elif curr_section is None:
                    raise NoSectionError(self.filename)

                else:
                    logger.info("File '%s' contains unexpected section '%s'",
                                self.filename, curr_section)
