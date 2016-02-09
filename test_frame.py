import unittest

from frame import Atom
from frame import Bead
from frame import Frame


class AtomTest(unittest.TestCase):
    def test_atom_create(self):
        atom = Atom(name="Name", num=0, type="Type")
        self.assertEqual("Name", atom.name)
        self.assertEqual(0, atom.num)
        self.assertEqual("Type", atom.type)


class BeadTest(unittest.TestCase):
    def test_bead_create(self):
        bead = Bead(name="Name", num=0, type="Type")
        self.assertEqual("Name", bead.name)
        self.assertEqual(0, bead.num)
        self.assertEqual("Type", bead.type)


class FrameTest(unittest.TestCase):
    def test_frame_create(self):
        frame = Frame()

    def test_frame_add_atoms(self):
        atom = Atom(name="Name", num=0, type="Type")
        frame = Frame()
        frame.addAtom(atom)
        self.assertEqual(atom, frame.atoms[0])
        self.assertTrue(atom is frame.atoms[0])


if __name__ == '__main__':
    unittest.main()
