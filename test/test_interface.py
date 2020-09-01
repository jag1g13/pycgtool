import unittest

from pycgtool.interface import Progress


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
