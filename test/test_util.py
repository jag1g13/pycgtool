import unittest

import numpy as np

from pycgtool.util import stat_moments


class UtilTest(unittest.TestCase):
    def test_stat_moments(self):
        t1 = [3, 3, 3, 3, 3]
        t2 = [1, 2, 3, 4, 5]
        np.testing.assert_allclose(np.array([3, 0]), stat_moments(t1))
        np.testing.assert_allclose(np.array([3, 2]), stat_moments(t2))


if __name__ == '__main__':
    unittest.main()
