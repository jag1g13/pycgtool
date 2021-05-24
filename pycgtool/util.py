"""This module contains some general purpose utility functions used in PyCGTOOL."""

import filecmp
import functools
import itertools
import logging
import os
import pathlib
import random
import typing

import mdtraj
import numpy as np

logger = logging.getLogger(__name__)


def circular_mean(values: typing.Iterable[float]) -> float:
    """Calculate average of angles on a circle.

    See https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    """
    if not isinstance(values, np.ndarray):
        values = np.array(values)

    x = np.cos(values)
    y = np.sin(values)
    vec = np.array([x, y]).T

    x_av, y_av = np.nanmean(vec, axis=0)
    return np.arctan2(y_av, x_av)


def circular_variance(values: typing.Iterable[float]) -> float:
    """Calculate variance of angles on a circle.

    See https://en.wikipedia.org/wiki/Mean_of_circular_quantities
    """
    if not isinstance(values, np.ndarray):
        values = np.array(values)

    two_pi = 2 * np.pi
    mean = circular_mean(values)

    diff = values - mean
    diff -= np.rint(diff / two_pi) * two_pi
    return np.nanmean(np.square(diff))


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
                        if node2.strip("+") == pair2[
                                0] and "+" + pair2[1] not in chain:
                            append_if_not_in(
                                ret, spare + (node1, node2, "+" + pair2[1]))
                        elif node2.strip("+") == pair2[
                                1] and "+" + pair2[0] not in chain:
                            append_if_not_in(
                                ret, spare + (node1, node2, "+" + pair2[0]))
            except AttributeError:
                pass

            # Reverse and look for links at the start of the chain
            chain = chain[::-1]

    return ret


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
    logger.warning('Existing file %s backed up as %s', name, new_name)
    return new_name


def sliding(vals: typing.Iterable):
    """Yield three values in a sliding window along an iterable.

    :param vals: Iterable to iterate over
    :return: Generator of tuples
    """
    it = iter(vals)
    prev = None

    try:
        current = next(it)

    except StopIteration:
        raise ValueError('Cannot make sliding window over empty iterable')

    for nxt in it:
        yield (prev, current, nxt)
        prev = current
        current = nxt

    yield (prev, current, None)


def cmp_file_whitespace_float(ref_filename,
                              test_filename,
                              rtol=0.01,
                              verbose=False):
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

        return cmp_whitespace_float(ref_lines,
                                    test_lines,
                                    rtol=rtol,
                                    verbose=verbose)


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
    for i, (ref_line, test_line) in enumerate(
            itertools.zip_longest(ref_lines, test_lines)):
        # Shortcut trivial comparisons
        if ref_line is None or test_line is None:
            diff_lines.append((i, ref_line, test_line))
            continue

        if ref_line == test_line:
            continue

        ref_toks = ref_line.split()
        test_toks = test_line.split()

        # Check for float comparison of tokens on line
        for ref_tok, test_tok in itertools.zip_longest(
                map(number_or_string, ref_toks),
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


def any_starts_with(iterable: typing.Iterable, char: str) -> bool:
    """Return True if any string(s) in nested iterable starts with a given character.

    :param iterable iterable: Nested iterable of strings to check
    :param str char: Char to check each string
    :return bool: True if any string in nested iterable starts with char, else False
    """
    if isinstance(iterable, str):
        return iterable.startswith(char)

    recurse = functools.partial(any_starts_with, char=char)
    return any(map(recurse, iterable))


def load_optional_topology(
    trajfile: typing.Union[str, pathlib.Path],
    topfile: typing.Optional[typing.Union[str, pathlib.Path]] = None
) -> mdtraj.Trajectory:
    """Load an MDTraj trajectory with or without a separate topology file.

    :param trajfile: Trajectory file
    :param topfile: Corresponding topology file
    :return: MDTraj trajectory
    """
    if topfile is None:
        return mdtraj.load(str(trajfile))

    return mdtraj.load(str(trajfile), top=str(topfile))


def compare_trajectories(
        *trajectory_files: typing.Union[str, pathlib.Path],
        topology_file: typing.Optional[typing.Union[str, pathlib.Path]] = None,
        rtol: float = 0.001) -> bool:
    """Compare multiple trajectory files for equality.

    :param trajectory_files: Paths of trajectory file to compare
    :param topology_file: Topology file path
    :param coords_only: Only compare coordinates from trajectory - e.g. not box size
    """
    ref_traj = load_optional_topology(trajectory_files[0], topology_file)

    for traj_file in trajectory_files[1:]:
        try:
            traj = load_optional_topology(traj_file, topology_file)

        except ValueError as exc:
            raise ValueError('Topology does not match') from exc

        # Conditions from Trajectory.__hash__
        if traj.topology != ref_traj.topology:
            raise ValueError('Topology does not match')
        if len(traj) != len(ref_traj):
            raise ValueError('Length does not match')
        if not np.allclose(traj.xyz, ref_traj.xyz, rtol=rtol):
            raise ValueError('Coordinates do not match')
        if not np.allclose(traj.time, ref_traj.time, rtol=rtol):
            raise ValueError('Time does not match')
        if not np.allclose(
                traj.unitcell_vectors, ref_traj.unitcell_vectors, rtol=rtol):
            raise ValueError('Unitcell does not match')

    return True
