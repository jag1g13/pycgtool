import numpy as np
np.seterr(all="raise")


def stat_moments(vals):
    """
    Return statistical (population) moments of data provided.

    :param vals: The data for which to calculate moments
    :return: Numpy array of moments - population mean and variance
    """
    res = np.zeros(2)
    try:
        for val in vals:
            res[0] += val
        mean = res[0] / len(vals)

        for val in vals:
            res[1] += pow(val - mean, 2)

        res /= len(vals)
        return res
    except FloatingPointError:
        return np.zeros(2)
