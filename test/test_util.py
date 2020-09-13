import unittest
import os
import logging
import tempfile
import shutil

import numpy as np

from pycgtool import util
from pycgtool.util import extend_graph_chain, transpose_and_sample
from pycgtool.util import dir_up, backup_file, sliding
from pycgtool.util import file_write_lines, cmp_whitespace_float
from pycgtool.util import circular_mean, circular_variance


class UtilTest(unittest.TestCase):
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

    def test_dir_up(self):
        path = os.path.realpath(__file__)
        self.assertEqual(path, dir_up(path, 0))
        self.assertEqual(os.path.dirname(path), dir_up(path))
        self.assertEqual(os.path.dirname(os.path.dirname(path)), dir_up(path, 2))

    def test_backup_file(self):
        try:
            os.remove("testfile")
            os.remove("#testfile.1#")
            os.remove("#testfile.2#")
        except OSError:
            pass

        logging.disable(logging.WARNING)
        open("testfile", "a").close()
        self.assertEqual("#testfile.1#", backup_file("testfile"))
        open("testfile", "a").close()
        self.assertTrue(os.path.exists("#testfile.1#"))
        self.assertEqual("#testfile.2#", backup_file("testfile"))
        open("testfile", "a").close()
        self.assertTrue(os.path.exists("#testfile.2#"))
        logging.disable(logging.NOTSET)

        os.remove("testfile")
        os.remove("#testfile.1#")
        os.remove("#testfile.2#")

    def test_sliding(self):
        l = [0, 1, 2, 3, 4]
        res = [(None, 0, 1), (0, 1, 2), (1, 2, 3), (2, 3, 4), (3, 4, None)]
        for res, pair in zip(res, sliding(l)):
            self.assertEqual(res, pair)

    def test_transpose_and_sample_no_sample(self):
        l = [(1, 2), (3, 4), (5, 6)]
        l_t = [(1, 3, 5), (2, 4, 6)]
        self.assertEqual(l_t, transpose_and_sample(l, None))

    def test_transpose_and_sample(self):
        l = [(1, 2), (3, 4), (5, 6)]
        l_t = [(1, 3, 5), (2, 4, 6)]

        l_t_test = transpose_and_sample(l, n=1)
        self.assertEqual(1, len(l_t_test))
        self.assertIn(l_t_test[0], l_t)

    def test_cmp_whitespace_text(self):
        ref = ["Hello World"]
        test = ["Hello World"]
        self.assertTrue(cmp_whitespace_float(ref, test))
        test = ["Hello PyCGTOOL"]
        self.assertFalse(cmp_whitespace_float(ref, test))

    def test_cmp_whitespace_float(self):
        ref =  ["8.3 1.00 -3 0"]
        test = ["8.3 1.00 -3 0.0"]
        self.assertTrue(cmp_whitespace_float(ref, test))
        test = ["8.3 1.01 -3 0"]
        self.assertFalse(cmp_whitespace_float(ref, test))
        self.assertTrue(cmp_whitespace_float(ref, test, rtol=0.1))

        test = ["8.3 1.00 -2.999 0"]
        self.assertTrue(cmp_whitespace_float(ref, test))
        test = ["8.3 1.00 -3.001 0"]
        self.assertTrue(cmp_whitespace_float(ref, test))

    def test_circular_mean(self):
        values = np.deg2rad(np.array([-175, 165], dtype=np.float32))
        test = np.deg2rad(175.)
        self.assertAlmostEqual(circular_mean(values),  test, delta=test / 500.)

        values = np.deg2rad(np.array([165, 145], dtype=np.float32))
        test = np.deg2rad(155.)
        self.assertAlmostEqual(circular_mean(values), test, delta=test / 500.)

    def test_circular_variance(self):
        values = np.deg2rad(np.array([165, 145], dtype=np.float32))
        test_var = np.var(values)
        self.assertAlmostEqual( circular_variance(values), test_var, delta=test_var / 500.)

        values = np.deg2rad(np.array([-175, 165], dtype=np.float32))
        self.assertAlmostEqual(circular_variance(values), test_var, delta=test_var / 500.)



# TODO test backing up
class UtilFileWriteLinesTest(unittest.TestCase):
    def setUp(self):
        self.test_dir = tempfile.mkdtemp()

    def tearDown(self):
        shutil.rmtree(self.test_dir)

    def test_file_write_lines_empty(self):
        filename = os.path.join(self.test_dir, "empty")
        file_write_lines(filename)
        with open(filename) as f:
            self.assertFalse(f.readlines())

    def test_file_write_lines_list(self):
        lines = ["a", "b c", "  ", "123 d", "", " %-+# "]

        filename = os.path.join(self.test_dir, "list")
        file_write_lines(filename, lines)
        with open(filename) as f:
            self.assertListEqual(lines, f.read().splitlines())

    def test_file_write_lines_append(self):
        lines = ["a", "b c", "d ", " ef"]

        filename = os.path.join(self.test_dir, "append")
        file_write_lines(filename, lines[:2])
        file_write_lines(filename, lines[2:], append=True)
        with open(filename) as f:
            self.assertListEqual(lines, f.read().splitlines())


class CompareTrajectoryTest(unittest.TestCase):
    def test_compare_trajectory_single(self):
        self.assertTrue(util.compare_trajectories(
            'test/data/sugar.gro',
            'test/data/sugar.gro',
            topology_file='test/data/sugar.gro'
        ))

    def test_compare_trajectory_single_false(self):
        with self.assertRaises(ValueError):
            util.compare_trajectories(
                'test/data/sugar.gro',
                'test/data/water.gro',
                topology_file='test/data/sugar.gro'
            )

    def test_compare_trajectory(self):
        self.assertTrue(util.compare_trajectories(
            'test/data/sugar.xtc',
            'test/data/sugar.xtc',
            topology_file='test/data/sugar.gro'
        ))

    def test_compare_trajectory_false(self):
        with self.assertRaises(ValueError):
            util.compare_trajectories(
                'test/data/sugar.xtc',
                'test/data/water.xtc',
                topology_file='test/data/sugar.gro'
            )


if __name__ == '__main__':
    unittest.main()
