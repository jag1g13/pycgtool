import unittest

from pycgtool.interface import Options, Progress


class UtilTest(unittest.TestCase):
    def test_progress_iterations(self):
        sum = 0
        for i in Progress(10, quiet=True):
            sum += i
        self.assertEqual(45, sum)

    @unittest.expectedFailure
    def test_progress_prewhile(self):
        sum = 0
        for i in Progress(20, prewhile=lambda: sum < 20, quiet=True):
           sum += i
        self.assertEqual(15, sum)

    def test_progress_postwhile(self):
        sum = 0
        for i in Progress(20, postwhile=lambda: sum < 20, quiet=True):
            sum += i
        self.assertEqual(21, sum)

    def test_options_init(self):
        opt = Options([("a", True), ("b", 10), ("c", "hello")])
        self.assertEqual(True, opt.a)
        self.assertEqual(10, opt.b)
        self.assertEqual("hello", opt.c)

    def test_options_set(self):
        opt = Options([("a", True), ("b", 10), ("c", "hello")])
        opt.set("a", False)
        opt.set("b", 11)
        self.assertEqual(False, opt.a)
        self.assertEqual(11, opt.b)
        self.assertEqual("hello", opt.c)
        opt.set("a", "yEs")
        opt.set("b", "12")
        self.assertEqual(True, opt.a)
        self.assertEqual(12, opt.b)
        opt.set("a", "f")
        self.assertEqual(False, opt.a)
        with self.assertRaises(ValueError):
            opt.set("a", "hello")
        with self.assertRaises(ValueError):
            opt.set("b", "hello")

    def test_options_set_by_num(self):
        opt = Options([("a", True), ("b", 10), ("c", "hello")])
        opt.set_by_num(0, False)
        opt.set_by_num(1, 11)
        self.assertEqual(False, opt.a)
        self.assertEqual(11, opt.b)
        self.assertEqual("hello", opt.c)
        opt.set_by_num(0, "yEs")
        opt.set_by_num(1, "12")
        self.assertEqual(True, opt.a)
        self.assertEqual(12, opt.b)
        opt.set_by_num(0, "f")
        self.assertEqual(False, opt.a)
        with self.assertRaises(ValueError):
            opt.set_by_num(0, "hello")
        with self.assertRaises(ValueError):
            opt.set_by_num(1, "hello")

    def test_options_toggle(self):
        opt = Options([("a", True), ("b", 10), ("c", "hello")])
        opt.toggle_boolean("a")
        self.assertEqual(False, opt.a)
        opt.toggle_boolean_by_num(0)
        self.assertEqual(True, opt.a)
