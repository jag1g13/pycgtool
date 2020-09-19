"""Module containing classes used to parse custom CFG file format.

Format is based upon GROMACS .itp files but does not support nesting of sections.
"""

import collections
import contextlib
import pathlib
import typing


class DuplicateSectionError(KeyError):
    """Exception used to indicate that a section has appeared twice in a file."""
    def __init__(self, section, filename):
        msg = f"Section '{section}' appears twice in file '{filename}'."
        super().__init__(msg)


class NoSectionError(KeyError):
    """Exception used to indicate that a file contains no sections."""
    def __init__(self, filename):
        msg = f"File '{filename}' contains no '[]' section headers."
        super().__init__(msg)


class CFG(collections.OrderedDict, contextlib.AbstractContextManager):
    """Class representing a CFG file.

    Contains a dictionary of Sections.
    """
    def __init__(self,
                 filename: typing.Optional[typing.Union[str,
                                                        pathlib.Path]] = None):
        """Parse a config file and extract Sections.

        :param filename: Name of file to read
        :return: Instance of CFG
        """
        super().__init__()
        self.filename = None

        if filename is not None:
            self.filename = pathlib.Path(filename)
            self._read_file()

    def _read_line(self, line: str) -> typing.Optional[str]:
        # Strip comments
        line = line.split(';')[0].strip()
        if not line:
            return None

        # Handle include directive
        if line.startswith('#include'):
            include_file = line.split()[1].strip('"')
            other = type(self)(self.filename.parent.joinpath(include_file))
            self.update(other)
            return None

        return line

    def _read_file(self) -> None:
        with open(self.filename) as f:
            curr_section = None

            for line in f:
                line = self._read_line(line)
                if line is None:
                    continue

                if line.startswith("["):
                    curr_section = line.strip("[ ]")
                    if curr_section in self:
                        raise DuplicateSectionError(curr_section,
                                                    self.filename)
                    self[curr_section] = []
                    continue

                toks = tuple(line.split())
                try:
                    self[curr_section].append(toks)
                except KeyError as e:
                    raise NoSectionError(self.filename) from e

    def __exit__(self, exc_type, exc_value, traceback):
        pass
