import unittest

from frame import Atom
from frame import Bead
from frame import Frame


class AtomTest(unittest.TestCase):
    def test_atom_create(self):
        atom = Atom(name="Name", num=0, type="Type")

class BeadTest(unittest.TestCase):
    def test_bead_create(self):
        bead = Bead(name="Name", num=0, type="Type")

class FrameTest(unittest.TestCase):
    def test_frame_create(self):
        frame = Frame()


if __name__ == '__main__':
    unittest.main()
