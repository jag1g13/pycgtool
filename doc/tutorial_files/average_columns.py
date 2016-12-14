#!/usr/bin/env python3
"""
Calculate mean and standard deviation of column data using Welford's one pass algorithm
"""

import sys
import numpy as np

np.set_printoptions(formatter={"float": lambda x: "{0:10.3f}".format(x)})

def average_columns(filename):
    with open(filename) as f:
        mean, sdev = None, None
        nlines = 0

        for line in f:
            if not line.strip():
                continue
            try:
                cols = np.array(list(map(float, line.split())))
            except ValueError:
                continue
            if mean is None:
                mean = np.zeros(len(cols))
                sdev = np.zeros(len(cols))
            nlines += 1
            delta = cols - mean
            mean += delta / nlines
            sdev += delta * (cols - mean)


        sdev = np.sqrt(sdev / (nlines - 1))
        print("File: {0}".format(filename))
        print("Mean: {0}".format(mean))
        print("Sdev: {0}".format(sdev))
        print()
    return (mean, sdev)

def compare_files(filename1, filename2):
    means = (average_columns(filename1)[0], average_columns(filename2)[0])
    print("Diff: {0}".format(means[1] - means[0]))
    print("%Dif: {0}".format(100. * (means[1] - means[0]) / means[0]))

if __name__ == "__main__":
    if len(sys.argv) == 2:
        average_columns(sys.argv[1])
    elif len(sys.argv) == 3:
        compare_files(sys.argv[1], sys.argv[2])
    else:
        print("Give one filename to average columns, or two filenames to compare.")
