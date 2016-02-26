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

    # TODO check this is correct
    @unittest.expectedFailure
    def test_r_squared(self):
        ref = [(i, i) for i in range(4)]
        fit = ref
        self.assertEqual(1, r_squared(ref, fit))
        fit = [(1.5, 1.5) for _ in range(4)]
        self.assertEqual(0, r_squared(ref, fit))
        fit = [(i, i) for i in range(3, -1, -1)]
        self.assertEqual(-1, r_squared(ref, fit))


if __name__ == '__main__':
    unittest.main()
