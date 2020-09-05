"""
Module containing classes used to parse ITP/TOP file format.
"""

import os
import collections

from .cfg import NoSectionError, DuplicateSectionError

class ITP(collections.OrderedDict):
    """
    Class representing an .itp file

    Contains a dictionary for every molecule definition
    """
    __slots__ = ["filename"]

    def __init__(self, filename=None):
        super(ITP, self).__init__()
        self.filename = filename
        mol_sections = ["atoms", "bonds", "angles", "dihedrals"]

        with open(self.filename) as f:
            curr_section = None
            curr_mol = None
            lines = f.readlines()
            for i, line in enumerate(lines):
                line = line.strip()
                if not line or line.startswith(";"):
                    continue

                elif line.startswith("#"):
                    if line.startswith("#include"):
                        print(line.split()[1].strip('"'), line)
                        itp2 = ITP(os.path.join(os.path.dirname(self.filename),
                                                line.split()[1].strip('"')))
                        self.update(itp2)
                    continue

                elif line.startswith("["):
                    curr_section = line.split(';')[0].strip("[ ]")
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

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass