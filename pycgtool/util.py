"""
This module contains some general purpose utility functions used in PyCGTOOL.
"""

from collections import namedtuple
import filecmp
import functools
import itertools
import logging
import math
import os
import random
import re
import typing

import mdtraj
import numpy as np

logger = logging.getLogger(__name__)


class NumbaDummy:
    """
    Dummy Numba module
    """
    def __getattr__(self, item):
        if item == "jit":
            return NumbaDummy.jit
        return self

    def __getitem__(self, item):
        return self

    def __call__(self, *args, **kwargs):
        return None

    @staticmethod
    def jit(*args, **kwargs):
        """
        Dummy version of numba.jit decorator, does nothing
        """
        if len(args) == 1 and callable(args[0]):
            return args[0]
        else:
            def wrap(func):
                return func
            return wrap

try:
    import numba
except ImportError:
    # If numba isn't installed, create dummy so we don't get NameErrors
    numba = NumbaDummy()


@numba.jit(numba.float32[3](numba.float32[3], numba.float32[3]))
def vector_cross(u, v):
    """
    Return vector cross product of two 3d vectors as numpy array.

    :param u: First 3d vector
    :param v: Second 3d vector
    :return: Cross product of two vectors as numpy.array
    """
    res = np.empty_like(u)
    res[0] = u[1] * v[2] - u[2] * v[1]
    res[1] = u[2] * v[0] - u[0] * v[2]
    res[2] = u[0] * v[1] - u[1] * v[0]
    return res

def circular_mean(values):
    """
    Return average of angles on a cirle

    See https://en.wikipedia.org/wiki/Mean_of_circular_quantities

    :param values: array values to average
    :return: average
    """
    if not isinstance(values, np.ndarray):
        values = np.array(values)

    x = np.cos(values)
    y = np.sin(values)
    vec = np.array([x, y]).T

    x_av, y_av = np.nanmean(vec, axis=0)
    return np.arctan2(y_av, x_av)

def circular_variance(values):
    """
    Return variance of angles on a cirle

    :param values: array values
    :return: average
    """
    if not isinstance(values, np.ndarray):
        values = np.array(values)

    two_pi = 2 * np.pi
    mean = circular_mean(values)

    diff = values - mean
    diff -= np.rint(diff / two_pi ) * two_pi
    return np.nanmean(np.square(diff))


@numba.jit(numba.float32(numba.float32[3], numba.float32[3]))
def vector_dot(u, v):
    """
    Return vector dot product of two 3d vectors.

    :param u: First 3d vector
    :param v: Second 3d vector
    :return: Dot product of two vectors
    """
    return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]


@numba.jit(numba.float32(numba.float32[3]))
def vector_len(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])


@numba.jit(numba.float32(numba.float32[3], numba.float32[3]))
def vector_angle(a, b):
    """
    Calculate the angle between two vectors.

    :param a: First vector
    :param b: Second vector
    :return: Angle in radians
    """
    mag = vector_len(a) * vector_len(b)
    if mag == 0:
        raise ZeroDivisionError("One or more bonds in angle calculation has length zero")
    dot = vector_dot(a, b) / mag
    # Previously checked to within 1%, but this prevented numba optimisation
    ang = math.acos(max(-1, min(1, dot)))
    return ang


@numba.jit
def vector_angle_signed(a, b, ref=np.array([0., 0., 1.])):
    """
    Calculate the signed angle between two vectors.

    :param a: First vector
    :param b: Second vector
    :param ref: Optional reference vector, will use global z-axis if not provided
    :return: Signed angle in radians
    """
    ang = vector_angle(a, b)
    signum = math.copysign(1, vector_dot(vector_cross(a, b), ref))
    return ang * signum


def tuple_equivalent(tuple1, tuple2):
    """
    Check if two node tuples are equivalent. Assumes undirected edges.

    :param tuple1: First tuple to compare
    :param tuple2: Second tuple to compare
    :return: True if tuples are equivalent, else False
    """
    if tuple1 == tuple2:
        return True
    elif tuple1[::-1] == tuple2:
        return True
    else:
        return False


@numba.jit
def dist_with_pbc(pos1, pos2, box):
    """
    Calculate the distance between two points accounting for periodicity.

    :param pos1: 3d vector position 1
    :param pos2: 3d vector position 2
    :param box: Cubic box vectors
    :return: Vector between two points
    """
    d = pos2 - pos1
    if box[0] * box[1] * box[2]:
        d -= box * np.rint(d / box)
    return d


def extend_graph_chain(extend, pairs):
    """
    Take list of tuples representing chained links in an undirected graph and extend the chain length.

    :param extend: List of link tuples to extend
    :param pairs: Graph edges as list of tuples
    :return: List of link tuples for chain length one greater than input
    """
    ret = []

    def append_if_not_in(lst, item):
        if item not in lst and item[::-1] not in lst:
            lst.append(item)

    for chain in extend:
        for _ in range(2):
            node1, node2 = chain[-2:]
            spare = chain[:-2]

            for pair2 in pairs:
                if node2 == pair2[0] and pair2[1] not in chain:
                    append_if_not_in(ret, spare + (node1, node2, pair2[1]))
                elif node2 == pair2[1] and pair2[0] not in chain:
                    append_if_not_in(ret, spare + (node1, node2, pair2[0]))

            try:
                # Support GROMACS RTP + to link to next residue
                if node2.startswith("+"):
                    for pair2 in pairs:
                        if node2.strip("+") == pair2[0] and "+" + pair2[1] not in chain:
                            append_if_not_in(ret, spare + (node1, node2, "+" + pair2[1]))
                        elif node2.strip("+") == pair2[1] and "+" + pair2[0] not in chain:
                            append_if_not_in(ret, spare + (node1, node2, "+" + pair2[0]))
            except AttributeError:
                pass

            # Reverse and look for links at the start of the chain
            chain = chain[::-1]

    return ret


def stat_moments(vals, ignore_nan=True):
    """
    Return statistical (population) moments of data provided.

    :param vals: The data for which to calculate moments
    :param ignore_nan: Whether to exclude np.nan and infinity from calculation
    :return: Tuple of moments - population mean and variance, both zero if input is empty list
    """
    if ignore_nan:
        vals_tmp = [val for val in vals if np.isfinite(val)]
    else:
        vals_tmp = vals

    try:
        return np.mean(vals_tmp), np.var(vals_tmp)
    except FloatingPointError:
        return 0., 0.


def transpose_and_sample(sequence, n=None):
    """
    Transpose a sequence of lists and sample to provide target number of rows.

    :param sequence: 2d sequence object to transpose
    :param n: Number of samples to take
    """
    rows = list(zip(*sequence))

    if n is not None and len(rows) > n:
        rows = random.sample(rows, n)

    return rows


def dir_up(name, n=1):
    """
    Return the directory path n levels above a specified file/directory.

    :param name: Name of file/directory to start from
    :param n: Number of directory levels to go up
    :return: Directory n directories above name
    """
    for _ in range(n):
        name = os.path.dirname(name)
    return name


def backup_file(name):
    """
    Backup a file using the GROMACS backup naming scheme.
    name -> #name.x#

    :param name: File to backup
    :return: New name of file after backup
    """
    if not os.path.exists(name):
        return None

    directory, basename = os.path.split(name)
    for i in itertools.count(1):
        new_base = "#{0}.{1}#".format(basename, i)
        new_name = os.path.join(directory, new_base)
        if not os.path.exists(new_name):
            break

    os.rename(name, new_name)
    logger.warning("Existing file {0} backed up as {1}".format(name, new_name))
    return new_name


def sliding(vals):
    """
    Yield three values in a sliding window along an iterable.

    :param vals: Iterable to iterate over
    :return: Generator of tuples
    """
    it = iter(vals)
    prev = None
    current = next(it)
    for nxt in it:
        yield (prev, current, nxt)
        prev = current
        current = nxt
    yield (prev, current, None)


def r_squared(ref, fit):
    """
    Calculate residual R squared of fitted data against reference data by y values.
    :param ref: Reference y values
    :param fit: Fitted y values
    :return: R squared
    """
    y_mean = sum(ref) / len(ref)
    ss_res, ss_tot = 0, 0
    for refpy, fitpy in zip(ref, fit):
        ss_res += (fitpy - refpy)**2
        ss_tot += (refpy - y_mean)**2
    try:
        return 1 - (ss_res / ss_tot)
    except ZeroDivisionError:
        return 0


def cmp_file_whitespace_float(ref_filename, test_filename, rtol=0.01, verbose=False):
    """
    Compare two files ignoring spacing on a line and using a tolerance on floats

    :param ref_filename: Name of reference file
    :param test_filename: Name of test file
    :param float rtol: Acceptable relative error on floating point numbers
    :param bool verbose: Print failing lines
    :return: True if files are the same, else False
    """
    if filecmp.cmp(ref_filename, test_filename):
        return True

    with open(ref_filename) as ref, open(test_filename) as test:
        ref_lines = ref.readlines()
        test_lines = test.readlines()

        return cmp_whitespace_float(ref_lines, test_lines, rtol=rtol, verbose=verbose)


def number_or_string(string: str) -> typing.Union[float, int, str]:
    """Convert string into an int or float if possible."""
    try:
        as_float = float(string)
        as_int = int(as_float)
        if as_int == as_float:
            return as_int

        return as_float

    except ValueError:
        return string


def cmp_whitespace_float(ref_lines, test_lines, rtol=0.01, verbose=False):
    """
    Compare two iterables of lines ignoring spacing on a line and using a tolerance on floats

    :param ref_lines: Iterable of reference lines
    :param test_lines: Iterable of test lines
    :param float rtol: Acceptable relative error on floating point numbers
    :param bool verbose: Print failing lines
    :return: True if all lines are the same, else False
    """
    diff_lines = []
    for i, (ref_line, test_line) in enumerate(itertools.zip_longest(ref_lines, test_lines)):
        # Shortcut trivial comparisons
        if ref_line is None or test_line is None:
            diff_lines.append((i, ref_line, test_line))
            continue

        if ref_line == test_line:
            continue

        ref_toks = ref_line.split()
        test_toks = test_line.split()

        if len(ref_toks) != len(test_toks):
            diff_lines.append((i, ref_line, test_line))

        # Check for float comparison
        for ref_tok, test_tok in zip(map(number_or_string, ref_toks),
                                     map(number_or_string, test_toks)):
            if ref_tok != test_tok:
                try:
                    if abs(ref_tok - test_tok) > abs(ref_tok) * rtol:
                        diff_lines.append((i, ref_line, test_line))
                except TypeError:
                    diff_lines.append((i, ref_line, test_line))

    if verbose and diff_lines:
        print("Lines fail comparison:")
        for i, ref_line, test_line in diff_lines:
            print("Line {}".format(i))
            print("Ref:  {0}".format(ref_line))
            print("Test: {0}".format(test_line))

    return len(diff_lines) == 0


def once_wrapper(func):
    """
    Wrap a function such that it runs only once, subsequent calls are ignored.

    :param func: Function to wrap
    :return: Wrapped function which will only run once.
    """
    called = False

    def wrap(*args, **kwargs):
        nonlocal called
        if not called:
            called = True
            return func(*args, **kwargs)

    return wrap


class SimpleEnum(object):
    class Enum(object):
        def __iter__(self):
            return iter(self.members)

        def __contains__(self, key):
            return key in self.members

        def __len__(self):
            return len(self.members)

        def __getitem__(self, key):
            return getattr(self, key)

        def as_dict(self):
            return {key: getattr(self, key).value for key in self.members}

    class EnumItem(object):
        def __init__(self, enum_name, key, value=None):
            self.enum_name = enum_name
            self.key = key
            if value is None:
                self.value = key
                self._default_value = True
            else:
                self.value = value
                self._default_value = False

        def __repr__(self):
            if self._default_value:
                return "<{0}: {1}>".format(self.enum_name, self.key)
            else:
                return "<{0}: {1} - {2}>".format(self.enum_name, self.key, self.value)

        def __hash__(self):
            return hash(repr(self))

        def __eq__(self, other):
            """
            Compare enums by key rather than by value.

            :param other: RHS enum
            :return: Whether EnumValues are equal by key
            """
            if self.enum_name != other.enum_name:
                raise TypeError("Cannot compare values from different enums")
            return self.key == other.key

        def compare_value(self, other):
            """
            Compare enums by value rather than by key.

            :param other: RHS enum
            :return: Whether EnumValues are equal by value
            """
            if self.enum_name != other.enum_name:
                raise TypeError("Cannot compare values from different enums")
            return self.value == other.value

    @classmethod
    def enum(cls, name, keys=None, values=None):
        def returner(val):
            return lambda _: val
        enum_cls = type(name, (cls.Enum,), {})

        if keys is None:
            keys = []
        if values is None:
            for key in keys:
                prop = property(returner(cls.EnumItem(name, key)))
                setattr(enum_cls, key, prop)
        else:
            for key, value in zip(keys, values):
                prop = property(returner(cls.EnumItem(name, key, value)))
                setattr(enum_cls, key, prop)

        setattr(enum_cls, "members", property(returner(set(keys))))
        return enum_cls()

    @classmethod
    def enum_from_dict(cls, name, key_val_dict):
        return cls.enum(name, key_val_dict.keys(), key_val_dict.values())


class FixedFormatUnpacker(object):
    """
    Unpack strings printed in fixed format.
    """
    FormatItem = namedtuple("FormatItem", ["type", "width"])
    FormatStyle = SimpleEnum.enum("FormatStyle", ["C", "Fortran"])

    def __init__(self, format_string, format_style=FormatStyle.C):
        """
        Does not support format strings containing spaces
        """
        self._format_string = format_string
        format_parsers = {self.FormatStyle.C:       self._parse_c_format,
                          self.FormatStyle.Fortran: self._parse_fortran_format}
        clean_format_string = self._clean_format_string(format_string)
        self._format_items = format_parsers[format_style](clean_format_string)

    @staticmethod
    def _clean_format_string(format_string):
        cleaned = format_string.lower()
        cleaned = cleaned.strip("-+#/")
        return cleaned

    @classmethod
    def _parse_c_format(cls, format_string):
        """
        Parse a C format specifier string.

        :param format_string: C format string
        :return: List of FormatItems
        """
        types = {"d": int, "s": str, "f": float}
        p = re.compile(r'%-?(?P<width>[0-9]+)[.0-9]*(?P<type>[dsfDSF])')
        items = [cls.FormatItem(types[match.group("type")], int(match.group("width"))) for match in p.finditer(format_string)]
        return items

    @classmethod
    def _parse_fortran_format(cls, format_string):
        """
        Parse a Fortran format specifier string.

        :param format_string: Fortran format string
        :return: List of FormatItems
        """
        items = []
        types = {"i": int, "a": str, "f": float, "x": None}
        p = re.compile(r'\s*(?P<repeat>[0-9]*)(?P<type>[ixafIXAF])(?P<width>[0-9]*)\s*')
        for match in p.finditer(format_string):
            repeat = int(match.group("repeat")) if match.group("repeat") else 1
            width = int(match.group("width")) if match.group("width") else 1
            for _ in range(repeat):
                items.append(cls.FormatItem(types[match.group("type")], width))
        return items

    def unpack(self, string):
        items = []
        start = 0
        for format_item in self._format_items:
            string_part = string[start:start+format_item.width]
            start += format_item.width
            if format_item.type is not None:
                items.append(format_item.type(string_part.strip()))
        return items


def tqdm_dummy(iterable, **kwargs):
    return iterable


def file_write_lines(filename, lines=None, backup=True, append=False):
    """
    Open a file and write lines to it.

    :param str filename: Name of file to open
    :param iterable[str] lines: Iterable of lines to write
    :param bool backup: Should the file be backed up if it exists?  Disabled if appending
    :param bool append: Should lines be appended to an existing file?
    """
    if backup and not append:
        backup_file(filename)

    mode = "a" if append else "w"
    with open(filename, mode) as f:
        if lines is not None:
            for line in lines:
                print(line, file=f)


def any_starts_with(iterable, char):
    """
    Return True if any string(s) in nested iterable starts with 'char'.

    :param iterable iterable: Nested iterable of strings to check
    :param str char: Char to check each string
    :return bool: True if any string in nested iterable starts with char, else False
    """
    recurse = functools.partial(any_starts_with, char=char)
    if type(iterable) is str:
        return iterable.startswith(char)
    else:
        return any(map(recurse, iterable))


def compare_trajectories(*trajectory_files: typing.Iterable[str],
                         topology_file: str,
                         coords_only: bool = False) -> bool:
    """Compare multiple trajectory files for equality.
    
    :param trajectory_files: Paths of trajectory file to compare
    :param topology_file: Topology file path
    :param coords_only: Only compare coordinates from trajectory - e.g. not box size
    """
    ref_traj = mdtraj.load(trajectory_files[0], top=topology_file)

    for traj_file in trajectory_files[1:]:
        try:
            traj = mdtraj.load(traj_file, top=topology_file)

        except ValueError:
            # Trajectory doesn't match topology
            return False

        if coords_only:
            if traj.xyz != ref_traj.xyz:
                return False

        else:
            if traj != ref_traj:
                return False

    return True
