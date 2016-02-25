"""
This module contains some general purpose utility functions used in PyCGTOOL.
"""

import numpy as np
np.seterr(all="raise")


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
