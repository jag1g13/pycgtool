"""
This module contains some general purpose utility functions used in PyCGTOOL.
"""

import numpy as np
np.seterr(all="raise")


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


def triplets_from_pairs(nodes, pairs):
    """
    Find all connected triplets from a list of connected pairs.
    Currently uses a naive algorithm, will be optimised if it turns out to be a performance bottleneck.

    :param nodes: List of node names
    :param pairs: List of node connectivity pairs
    :return: List of node connectivity triplets
    """
    triplets = []

    for node1 in nodes:
        for node2 in nodes:
            pair1 = (node1, node2)
            if pair1 not in pairs and pair1[::-1] not in pairs:
                continue

            for node3 in nodes:
                if node3 == node1:
                    continue
                pair2 = (node2, node3)
                if pair2 not in pairs and pair2[::-1] not in pairs:
                    continue

                triplet = (node1, node2, node3)
                if triplet not in triplets and triplet[::-1] not in triplets:
                    triplets.append(triplet)

    return triplets


def quadruplets_from_pairs(nodes, pairs):
    """
    Find all connected quadruplets from a list of connected pairs.
    Currently uses a naive algorithm, will be optimised if it turns out to be a performance bottleneck.

    :param nodes: List of node names
    :param pairs: List of node connectivity pairs
    :return: List of node connectivity quadruplets
    """
    quadruplets = []

    for node1 in nodes:
        for node2 in nodes:
            pair1 = (node1, node2)
            if pair1 not in pairs and pair1[::-1] not in pairs:
                continue

            for node3 in nodes:
                if node3 == node1:
                    continue
                pair2 = (node2, node3)
                if pair2 not in pairs and pair2[::-1] not in pairs:
                    continue

                for node4 in nodes:
                    if node4 == node2 or node4 == node1:
                        continue
                    pair3 = (node3, node4)
                    if pair3 not in pairs and pair3[::-1] not in pairs:
                        continue

                    quadruplet = (node1, node2, node3, node4)
                    if quadruplet not in quadruplets and quadruplet[::-1] not in quadruplets:
                        quadruplets.append(quadruplet)

    return quadruplets


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
    gaussian_single = lambda x: (amplitude / (sdev * np.sqrt(2*np.pi))) * np.exp((-(x-mean)**2) / (2*sdev*sdev))
    return map(gaussian_single, xs)
