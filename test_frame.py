import unittest

from frame import *


class MappingTest(unittest.TestCase):
    def test_mapping_create(self):
        mappings = mappings_from_file("data/water.map")
        self.assertEqual(1, len(mappings))
        self.assertEqual(3, len(mappings["SOL"].atoms))
        self.assertEqual("OW", len(mappings["SOL"].atoms[0].name))


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


class ResidueTest(unittest.TestCase):
    def test_residue_create(self):
        residue = Residue(name="Resname")
        self.assertEqual("Resname", residue.name)

    def test_residue_add_atoms(self):
        atom = Atom(name="Name", num=0, type="Type")
        residue = Residue()
        residue.add_atom(atom)
        self.assertEqual(atom, residue.atoms[0])
        self.assertTrue(atom is residue.atoms[0])


class FrameTest(unittest.TestCase):
    def test_frame_create(self):
        frame = Frame()

    def test_frame_add_residue(self):
        residue = Residue()
        frame = Frame()
        frame.add_residue(residue)
        self.assertEqual(residue, frame.residues[0])
        self.assertTrue(residue is frame.residues[0])

    def test_frame_read_gro(self):
        frame = Frame(gro="data/water.gro")
        self.assertEqual(4074, len(frame.residues))
        self.assertEqual("SOL", frame.residues[0].name)
        self.assertEqual(3, len(frame.residues[0].atoms))
        self.assertEqual("OW", frame.residues[0].atoms[0].name)


if __name__ == '__main__':
    unittest.main()
