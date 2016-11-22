"""
This module contains some general purpose utility functions used in PyCGTOOL.
"""

import os
import itertools
import random
import math
import filecmp
import logging

from collections import namedtuple

import numpy as np
np.seterr(all="raise")

logger = logging.getLogger(__name__)


class NumbaDummy:
    """
    Dummy Numba module
    """
    def __getattr__(self, item):
        if item == "jit":
            return self.jit
        return self

    def __getitem__(self, item):
        return self

    def __call__(self, *args, **kwargs):
        return None

    def jit(self, *args, **kwargs):
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
    if box[0] * box[1] * box[2] != 0:
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


def cmp_whitespace_float(ref_filename, test_filename, float_rel_error=0.01):
    """
    Compare two files ignoring spacing on a line and using a tolerance on floats

    :param ref_filename: Name of reference file
    :param test_filename: Name of test file
    :param float_rel_error: Acceptable relative error on floating point numbers
    :return: True if files are the same, else False
    """
    def float_or_string(string):
        try:
            int(string)
            return string
        except ValueError:
            try:
                return float(string)
            except ValueError:
                return string

    if filecmp.cmp(ref_filename, test_filename):
        return True
    with open(ref_filename) as ref_file, open(test_filename) as test_file:
        for ref_line, test_line in itertools.zip_longest(ref_file, test_file):
            if ref_line is None or test_line is None:
                return False
            if ref_line == test_line:
                continue

            ref_toks = ref_line.split()
            test_toks = test_line.split()
            if len(ref_toks) != len(test_toks):
                return False
            for ref_tok, test_tok in zip(map(float_or_string, ref_toks),
                                         map(float_or_string, test_toks)):
                if ref_tok != test_tok:
                    if float == type(ref_tok) and float == type(test_tok):
                        if abs(ref_tok - test_tok) > ref_tok * float_rel_error:
                            return False
                    else:
                        return False
    return True


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


# TODO testing
class FixedFormatUnpacker(object):
    """
    Unpack strings printed in fixed format.
    """
    types = {"d": int, "s": str, "f": float}
    FormatItem = namedtuple("FormatItem", ["type", "width"])

    def __init__(self, format_string):
        """
        Does not support format strings containing spaces
        """
        self._format_string = format_string
        self._format_items = []
        for item in format_string.split("%")[1:]:
            item_width = int(item[:-1])
            item_type = self.types[item[-1]]
            self._format_items.append(self.FormatItem(item_type, item_width))
    
    def unpack(self, string):
        items = []
        start = 0
        for format_item in self._format_items:
            string_part = string[start:start+format_item.width]
            start += format_item.width
            items.append(format_item.type(string_part.strip()))
        return items
