"""
Module containing classes used to parse custom CFG file format.

Format is based upon GROMACS .itp files but does not support nesting of sections.
"""
import os

from collections import OrderedDict


class Section:
    """
    Class representing a single section of the config file.
    """
    __slots__ = ["name", "_lines"]

    def __init__(self, name=None):
        """
        Create a section and storage for the lines it contains.

        :param name: Name of section
        """
        self.name = name
        self._lines = []

    def __len__(self):
        return len(self._lines)

    def __iter__(self):
        return iter(self._lines)

    def __getitem__(self, item):
        return self._lines[item]

    def add_line(self, line):
        self._lines.append(line)


class DuplicateSectionError(Exception):
    """
    Exception used to indicate that a section has appeared twice in a file.
    """
    def __repr__(self):
        return "Section {0} appears twice in file {1}.".format(*self.args)


class CFG:
    """
    Class representing a CFG file.

    Contains a dictionary of Sections.
    """
    __slots__ = ["filename", "_sections"]

    def __init__(self, filename=None, allow_duplicate=False):
        """
        Parse a config file and extract Sections.

        :param filename: Name of file to read
        :param allow_duplicate: Allow sections to appear more than once in a file
        :return: Instance of CFG
        """
        self.filename = filename
        self._sections = OrderedDict()

        with open(self.filename) as f:
            curr_section = None
            for line in f:
                line = line.strip()
                if not line or line.startswith(";"):
                    continue

                elif line.startswith("#include"):
                    cfg2 = CFG(os.path.join(os.path.dirname(self.filename),
                                            line.split()[1]))
                    self._sections.update(cfg2._sections)
                    continue

                elif line.startswith("["):
                    curr_section = line.strip("[ ]")
                    if curr_section in self._sections and not allow_duplicate:
                        raise DuplicateSectionError(curr_section, self.filename)
                    self._sections[curr_section] = Section(name=curr_section)
                    continue

                toks = tuple(line.split())
                self._sections[curr_section].add_line(toks)

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def __len__(self):
        return len(self._sections)

    def __iter__(self):
        return iter(self._sections.values())

    def __contains__(self, item):
        return item in self._sections

    def __getitem__(self, item):
        return self._sections[item]
