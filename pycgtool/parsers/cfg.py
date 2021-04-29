"""Module containing classes used to parse custom CFG file format.

Format is based upon GROMACS .itp files but does not support nesting of sections.
"""

import collections
import contextlib
import pathlib
import typing

PathLike = typing.Union[str, pathlib.Path]


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
    """Class representing a CFG file using a format similar to a GROMACS .itp file.

    Contains a dictionary of Sections.
    """
    def __init__(self, filepath: typing.Optional[PathLike] = None):
        """Parse a config file and extract Sections."""
        super().__init__()
        self.filepath = None

        if filepath is not None:
            self.filepath = pathlib.Path(filepath)
            self._read_file(self.filepath)

    def _read_line(self, line: str, filepath: pathlib.Path) -> str:
        # Strip comments
        line = line.split(';')[0].strip()

        # Handle include directive
        if line.startswith('#include'):
            include_file = line.split()[1].strip('"')
            other = type(self)(filepath.parent.joinpath(include_file))
            self.update(other)

            return ''  # Handle include then treat as empty line

        return line

    def _read_file(self, filepath: pathlib.Path) -> None:
        with open(filepath) as cfg_file:
            curr_section = None

            for line in cfg_file:
                line = self._read_line(line, filepath)

                if not line:
                    continue

                if line.startswith("["):
                    curr_section = line.strip("[ ]")

                    if curr_section in self:
                        raise DuplicateSectionError(curr_section, filepath)

                    self[curr_section] = []

                    continue

                toks = tuple(line.split())

                try:
                    self[curr_section].append(toks)

                except KeyError as exc:
                    raise NoSectionError(filepath) from exc

    def __exit__(self, exc_type, exc_value, traceback):
        pass
