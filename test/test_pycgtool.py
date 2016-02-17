import unittest
import subprocess
import os


class PycgtoolTest(unittest.TestCase):
    def test_run(self):
        path = os.path.dirname(os.path.dirname(__file__))
        self.assertEqual(0, subprocess.check_call([os.path.join(path, "pycgtool.py"), "-h"], stdout=subprocess.PIPE))


if __name__ == '__main__':
    unittest.main()
