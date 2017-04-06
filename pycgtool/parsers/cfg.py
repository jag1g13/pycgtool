"""
Module containing classes used to parse custom CFG file format.

Format is based upon GROMACS .itp files but does not support nesting of sections.
"""
import os

import collections


class DuplicateSectionError(KeyError):
    """
    Exception used to indicate that a section has appeared twice in a file.
    """
    def __init__(self, section, filename):
        msg = "Section '{0}' appears twice in file '{1}'.".format(section, filename)
        super(DuplicateSectionError, self).__init__(msg)


class NoSectionError(KeyError):
    """
    Exception used to indicate that a file contains no sections.
    """
    def __init__(self, filename):
        msg = "File '{0}' contains no '[]' section headers.".format(filename)
        super(NoSectionError, self).__init__(msg)


class CFG(collections.OrderedDict):
    """
    Class representing a CFG file.

    Contains a dictionary of Sections.
    """
    __slots__ = ["filename"]

    def __init__(self, filename=None):
        """
        Parse a config file and extract Sections.

        :param filename: Name of file to read
        :return: Instance of CFG
        """
        super(CFG, self).__init__()
        self.filename = filename

        with open(self.filename) as f:
            curr_section = None
            for line in f:
                line = line.strip()
                if not line or line.startswith(";"):
                    continue

                elif line.startswith("#include"):
                    cfg2 = CFG(os.path.join(os.path.dirname(self.filename),
                                            line.split()[1]))
                    self.update(cfg2)
                    continue

                elif line.startswith("["):
                    curr_section = line.strip("[ ]")
                    if curr_section in self:
                        raise DuplicateSectionError(curr_section, self.filename)
                    self[curr_section] = []
                    continue

                toks = tuple(line.split())
                try:
                    self[curr_section].append(toks)
                except KeyError as e:
                    raise NoSectionError(self.filename) from e

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass
