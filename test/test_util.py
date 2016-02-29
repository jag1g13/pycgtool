import unittest

import numpy as np

from pycgtool.util import stat_moments, sliding, r_squared


class UtilTest(unittest.TestCase):
    def test_stat_moments(self):
        t1 = [3, 3, 3, 3, 3]
        t2 = [1, 2, 3, 4, 5]
        np.testing.assert_allclose(np.array([3, 0]), stat_moments(t1))
        np.testing.assert_allclose(np.array([3, 2]), stat_moments(t2))

    def test_sliding(self):
        l = [0, 1, 2, 3, 4]
        res = [(None, 0, 1), (0, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, None)]
        for res, pair in zip(res, sliding(l)):
            self.assertEqual(res, pair)

    def test_r_squared(self):
        ref = [i for i in range(5)]
        fit = ref
        self.assertEqual(1, r_squared(ref, fit))
        fit = [2 for _ in range(5)]
        self.assertEqual(0, r_squared(ref, fit))
        fit = [i for i in range(1, 6)]
        self.assertEqual(0.5, r_squared(ref, fit))


if __name__ == '__main__':
    unittest.main()
