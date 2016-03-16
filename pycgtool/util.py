"""
This module contains some general purpose utility functions used in PyCGTOOL.
"""

import os
import itertools
import collections

import numpy as np
np.seterr(all="raise")


class Options:
    __slots__ = ['_dict']

    def __init__(self, default):
        self._dict = collections.OrderedDict()
        for key, val in default:
            try:
                val = val.lower()
            except AttributeError:
                pass

            self._dict[key.lower()] = (val, type(val))

    def __getattr__(self, attr):
        return self._dict[attr.lower()][0]

    def __repr__(self):
        res = "[" + ", ".join((str((key, val[0])) for key, val in self._dict.items())) + "]"
        return res

    def __iter__(self):
        return iter(((key, val[0]) for key, val in self._dict.items()))

    def __len__(self):
        return len(self._dict)

    def __getitem__(self, item):
        try:
            return self._dict[item]
        except KeyError:
            try:
                opt = list(self._dict.keys())[item]
                return self._dict[opt][0]
            except TypeError:
                raise TypeError("Must access Options using either a string or an integer")

    def set(self, opt, val):
        opt = opt.lower()
        try:
            val = val.lower()
        except AttributeError:
            pass
        _type = self._dict[opt][1]

        if _type is not type(val):
            if _type is bool:
                self._dict[opt] = (truthy(val), bool)
            else:
                self._dict[opt] = (_type(val), _type)
        else:
            self._dict[opt] = (val, _type)

    def set_by_num(self, opt_num, val):
        opt = list(self._dict.keys())[opt_num]
        self.set(opt, val)

    def toggle_boolean(self, opt):
        entry = self._dict[opt]
        if entry[1] is bool:
            self._dict[opt] = (not entry[0], entry[1])
        else:
            raise TypeError("Only boolean options can be toggled")

    def toggle_boolean_by_num(self, opt_num):
        opt = list(self._dict.keys())[opt_num]
        self.toggle_boolean(opt)


class Progress:
    def __init__(self, maxits, length=20, prewhile=None, postwhile=None):
        """
        Class to handle printing of a progress bar within loops.

        :param maxits: Expected number of iterations
        :param length: Length of progress bar in characters
        """
        self._maxits = maxits
        self._length = length
        self._prewhile = prewhile
        self._postwhile = postwhile
        self._its = 0

    def __enter__(self):
        return self

    def __exit__(self, type, value, traceback):
        pass

    def __iter__(self):
        return self

    def __next__(self):
        if self._postwhile is not None and self._its > 0 and not self._postwhile():
            self._stop()

        if self._prewhile is not None and not self._prewhile():
            self._stop()

        self._its += 1
        if self._its % 10 == 0:
            self._display()

        if self._its > self._maxits:
            self._stop()

        return self._its

    def _stop(self):
        self._display()
        print()
        raise StopIteration

    def _display(self):
        done = int(self._length * (self._its / self._maxits))
        left = self._length - done
        print("\r {0} [".format(self._its) + done * "#" + left * "-" + "] {0}".format(self._maxits), end="")


def truthy(string):
    truthy_strings = ("yes", "y", "on", "true", "t", "1")
    falsey_strings = ("no", "n", "off", "false", "f", "0")

    string = string.lower().strip()
    if string in truthy_strings:
        return True
    elif string in falsey_strings:
        return False
    else:
        raise ValueError("Value '{0}' could not be converted to boolean".format(string))


def cross(u, v):
    """
    Return vector cross product of two 3d vectors as numpy array.

    :param u: First 3d vector
    :param v: Second 3d vector
    :return: Cross product of two vectors as numpy.array
    """
    res = np.zeros(3)
    res[0] = u[1] * v[2] - u[2] * v[1]
    res[1] = u[2] * v[0] - u[0] * v[2]
    res[2] = u[0] * v[1] - u[1] * v[0]
    return res


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
    :param ignore_nan: Whether to exclude np.nan from calculation
    :return: Numpy array of moments - population mean and variance
    """
    if ignore_nan:
        vals_tmp = [val for val in vals if np.isfinite(val)]
    else:
        vals_tmp = vals

    res = np.zeros(2)
    try:
        for val in vals_tmp:
            res[0] += val
        mean = res[0] / len(vals_tmp)

        for val in vals_tmp:
            res[1] += pow(val - mean, 2)

        res /= len(vals_tmp)
        return res
    except FloatingPointError:
        return np.zeros(2)


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


def backup_file(name, verbose=False):
    """
    Backup a file using the GROMACS backup naming scheme.
    name -> #name.x#

    :param name: File to backup
    :param verbose: Print a message if file is backed up
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
    if verbose:
        print("Existing file {0} backed up as {1}".format(name, new_name))
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


def gaussian(xs, mean=0, sdev=1, amplitude=1):
    """
    Return y values from a Gaussian/normal distribution with provided mean, standard deviation and amplitude at x coordinates in vector x.

    :param xs: X values at which to calculate y values
    :param mean: Mean of Gaussian distribution
    :param sdev: Standard devaition of Gaussian distribution
    :param amplitude: Amplitude of Gaussian distribuion
    :return: Y values of Gaussian distribution at X values in x
    """
    def gaussian_single(x):
        prefactor = amplitude / (sdev * np.sqrt(2 * np.pi))
        top_bit = (x - mean) * (x - mean) / (2 * sdev * sdev)
        return prefactor * np.exp(-top_bit)
    return map(gaussian_single, xs)
