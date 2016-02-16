import numpy as np


def stat_moments(vals):
    """
    Return statistical moments of data provided.

    Note: results are population moments, not sample moments

    :param vals: The data for which to calculate moments
    :param moments: The number of moments required - default 2 (mean, var)
    :return: Numpy array of moments
    """
    res = np.zeros(2)
    for val in vals:
        res[0] += val
    mean = res[0] / len(vals)

    for val in vals:
        res[1] += pow(val - mean, 2)

    return res / len(vals)

