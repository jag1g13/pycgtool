import unittest

from pycgtool.interface import Options, Progress


class UtilTest(unittest.TestCase):
    def test_progress_iterations(self):
        sum = 0
        for i in Progress(10, quiet=True):
            sum += i
        self.assertEqual(45, sum)

    def test_progress_postwhile(self):
        sum = 0
        for i in Progress(20, dowhile=lambda: sum < 20, quiet=True):
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

    def test_options_toggle(self):
        opt = Options([("a", True), ("b", 10), ("c", "hello")])
        opt.toggle_boolean("a")
        self.assertEqual(False, opt.a)
