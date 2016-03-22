import unittest
import os

import numpy as np
import numpy.testing

from pycgtool.util import tuple_equivalent, extend_graph_chain, stat_moments
from pycgtool.util import dir_up, backup_file, sliding, r_squared, dist_with_pbc


class UtilTest(unittest.TestCase):
    def test_tuple_equivalent(self):
        t1 = (0, 1, 2)
        t2 = (0, 1, 2)
        self.assertTrue(tuple_equivalent(t1, t2))
        t2 = (2, 1, 0)
        self.assertTrue(tuple_equivalent(t1, t2))
        t2 = (2, 1, 3)
        self.assertFalse(tuple_equivalent(t1, t2))

    def test_dist_with_pbc(self):
        pos_a = np.array([1., 1., 1.])
        pos_b = np.array([9., 9., 9.])
        numpy.testing.assert_equal(np.array([8., 8., 8.]),
                                   dist_with_pbc(pos_a, pos_b, np.array([0., 0., 0.])))
        numpy.testing.assert_equal(np.array([8., 8., 8.]),
                                   dist_with_pbc(pos_a, pos_b, np.array([20., 20., 20.])))
        numpy.testing.assert_equal(np.array([-2., -2., -2.]),
                                   dist_with_pbc(pos_a, pos_b, np.array([10., 10., 10.])))

    def test_triplets_from_pairs(self):
        pairs = [(0, 1), (1, 2), (2, 3)]
        result = [(0, 1, 2), (1, 2, 3)]
        self.assertEqual(result, sorted(extend_graph_chain(pairs, pairs)))
        pairs = [(0, 1), (1, 2), (2, 3), (3, 0)]
        result = [(0, 1, 2), (1, 0, 3), (1, 2, 3), (2, 3, 0)]
        self.assertEqual(result, sorted(extend_graph_chain(pairs, pairs)))

    def test_triplets_from_pairs_multires(self):
        pairs = [("a", "b"), ("b", "c"), ("c", "d"), ("d", "+a")]
        result = [("a", "b", "c"), ("b", "c", "d"), ("c", "d", "+a"), ("d", "+a", "+b")]
        self.assertEqual(result, sorted(extend_graph_chain(pairs, pairs)))

    def test_quadruplets_from_pairs(self):
        pairs = [(0, 1), (1, 2), (2, 3)]
        result = [(0, 1, 2, 3)]
        triplets = extend_graph_chain(pairs, pairs)
        self.assertEqual(result, sorted(extend_graph_chain(triplets, pairs)))
        pairs = [(0, 1), (1, 2), (2, 3), (3, 0)]
        triplets = extend_graph_chain(pairs, pairs)
        result = [(0, 1, 2, 3), (1, 0, 3, 2), (1, 2, 3, 0), (2, 1, 0, 3)]
        self.assertEqual(result, sorted(extend_graph_chain(triplets, pairs)))

    def test_stat_moments(self):
        t1 = [3, 3, 3, 3, 3]
        t2 = [1, 2, 3, 4, 5]
        np.testing.assert_allclose(np.array([3, 0]), stat_moments(t1))
        np.testing.assert_allclose(np.array([3, 2]), stat_moments(t2))

    def test_dir_up(self):
        path = os.path.realpath(__file__)
        self.assertEqual(path, dir_up(path, 0))
        self.assertEqual(os.path.dirname(path), dir_up(path))
        self.assertEqual(os.path.dirname(os.path.dirname(path)), dir_up(path, 2))

    def test_backup_file(self):
        open("testfile", "a").close()
        self.assertEqual("#testfile.1#", backup_file("testfile"))
        open("testfile", "a").close()
        self.assertTrue(os.path.exists("#testfile.1#"))
        self.assertEqual("#testfile.2#", backup_file("testfile"))
        open("testfile", "a").close()
        self.assertTrue(os.path.exists("#testfile.2#"))
        os.remove("testfile")
        os.remove("#testfile.1#")
        os.remove("#testfile.2#")

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
