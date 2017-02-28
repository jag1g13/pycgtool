import unittest
import os
import logging

import numpy as np
import numpy.testing

from pycgtool.util import tuple_equivalent, extend_graph_chain, stat_moments, transpose_and_sample
from pycgtool.util import dir_up, backup_file, sliding, r_squared, dist_with_pbc
from pycgtool.util import SimpleEnum, FixedFormatUnpacker


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

    def test_r_squared(self):
        ref = [i for i in range(5)]
        fit = ref
        self.assertEqual(1, r_squared(ref, fit))
        fit = [2 for _ in range(5)]
        self.assertEqual(0, r_squared(ref, fit))
        fit = [i for i in range(1, 6)]
        self.assertEqual(0.5, r_squared(ref, fit))

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

    def test_simple_enum(self):
        enum = SimpleEnum.enum("enum", ["one", "two", "three"])
        self.assertTrue(enum.one == enum.one)
        self.assertTrue(enum.one.compare_value(enum.one))

        self.assertFalse(enum.two == enum.three)
        self.assertFalse(enum.two.compare_value(enum.three))

        with self.assertRaises(AttributeError):
            _ = enum.four
        with self.assertRaises(AttributeError):
            enum.one = 2

        enum2 = SimpleEnum.enum("enum2", ["one", "two", "three"])
        with self.assertRaises(TypeError):
            assert enum2.one == enum.one

        self.assertTrue("one" in enum)
        self.assertFalse("four" in enum)

    def test_simple_enum_values(self):
        enum = SimpleEnum.enum_from_dict("enum", {"one": 111,
                                                  "two": 111,
                                                  "three": 333})
        self.assertTrue(enum.one == enum.one)
        self.assertTrue(enum.one.compare_value(enum.one))

        self.assertFalse(enum.one == enum.two)
        self.assertTrue(enum.one.compare_value(enum.two))

        self.assertFalse(enum.two == enum.three)
        self.assertFalse(enum.two.compare_value(enum.three))

        with self.assertRaises(AttributeError):
            _ = enum.four
        with self.assertRaises(AttributeError):
            enum.one = 2

        enum2 = SimpleEnum.enum("enum2", ["one", "two", "three"])
        with self.assertRaises(TypeError):
            assert enum2.one == enum.one

        self.assertTrue("one" in enum)
        self.assertEqual(111, enum.one.value)

        self.assertFalse("four" in enum)

    def test_fixed_format_unpacker_c(self):
        unpacker = FixedFormatUnpacker("%-4d%5s%4.1f")
        toks = unpacker.unpack("1234hello12.3")
        self.assertEqual(3, len(toks))
        self.assertEqual(1234, toks[0])
        self.assertEqual("hello", toks[1])
        self.assertAlmostEqual(12.3, toks[2])

    def test_fixed_format_unpacker_fortran(self):
        unpacker = FixedFormatUnpacker("I4,A5,F4.1",
                                       FixedFormatUnpacker.FormatStyle.Fortran)
        toks = unpacker.unpack("1234hello12.3")
        self.assertEqual(3, len(toks))
        self.assertEqual(1234, toks[0])
        self.assertEqual("hello", toks[1])
        self.assertAlmostEqual(12.3, toks[2])

    def test_fixed_format_unpacker_fortran_space(self):
        unpacker = FixedFormatUnpacker("I4,X3,A5,X2,F4.1",
                                       FixedFormatUnpacker.FormatStyle.Fortran)
        toks = unpacker.unpack("1234 x hello x12.3")
        self.assertEqual(3, len(toks))
        self.assertEqual(1234, toks[0])
        self.assertEqual("hello", toks[1])
        self.assertAlmostEqual(12.3, toks[2])

    def test_fixed_format_unpacker_fortran_repeat(self):
        unpacker = FixedFormatUnpacker("2I2,X3,A5,2X,F4.1",
                                       FixedFormatUnpacker.FormatStyle.Fortran)
        toks = unpacker.unpack("1234 x hello x12.3")
        self.assertEqual(4, len(toks))
        self.assertEqual(12, toks[0])
        self.assertEqual(34, toks[1])
        self.assertEqual("hello", toks[2])
        self.assertAlmostEqual(12.3, toks[3])


if __name__ == '__main__':
    unittest.main()
