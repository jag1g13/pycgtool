import filecmp
import os
import pathlib
import unittest

import numpy as np

from pycgtool.mapping import BeadMap, Mapping, VirtualMap
from pycgtool.frame import Frame


class DummyOptions:
    map_center = "geom"
    virtual_map_center = "geom"


class BeadMapTest(unittest.TestCase):
    def test_beadmap_create(self):
        """Test that a BeadMap can be created."""
        bead = BeadMap("BEAD", 1, atoms=["ATOM1", "ATOM2"])

        self.assertEqual("BEAD", bead.name)
        self.assertEqual(1, bead.num)

        self.assertEqual(2, len(bead))
        self.assertEqual("ATOM1", bead[0])
        self.assertEqual("ATOM2", bead[1])

    def test_beadmap_guess_masses(self):
        """Test that masses can be assigned from atom names."""
        bead = BeadMap("BEAD", 1, atoms=["C1", "H1", "H2"])
        bead.guess_atom_masses()

        # Total mass of bead
        self.assertEqual(14, int(bead.mass))

        # Contribution of each atom to bead mass
        np.testing.assert_allclose(
            np.array([12, 1, 1]) / 14, bead.weights_dict["mass"], rtol=0.01
        )  # Account for 1% error in approximate masses


class MappingTest(unittest.TestCase):
    base_dir = pathlib.Path(__file__).absolute().parent
    data_dir = base_dir.joinpath("data")

    def test_mapping_create(self):
        """Test that a mapping can be correctly read from a file."""
        mapping = Mapping(self.data_dir.joinpath("water.map"), DummyOptions)
        self.assertEqual(2, len(mapping))  # SOL and HOH
        self.assertTrue("SOL" in mapping)

        sol_map = mapping["SOL"]

        self.assertEqual(1, len(sol_map))
        self.assertEqual(3, len(sol_map[0].atoms))
        self.assertEqual("OW", sol_map[0].atoms[0])
        self.assertEqual("HW1", sol_map[0].atoms[1])
        self.assertEqual("HW2", sol_map[0].atoms[2])

    def test_mapping_rename(self):
        """Test that an alternate mapping is created using MDTraj conventions."""
        mapping = Mapping(self.data_dir.joinpath("water.map"), DummyOptions)
        self.assertEqual(2, len(mapping))  # SOL and HOH
        self.assertTrue("HOH" in mapping)

        sol_map = mapping["HOH"]

        self.assertEqual(1, len(sol_map))
        self.assertEqual(3, len(sol_map[0].atoms))
        self.assertEqual("O", sol_map[0].atoms[0])
        self.assertEqual("H1", sol_map[0].atoms[1])
        self.assertEqual("H2", sol_map[0].atoms[2])

    def test_virtual_mapping_create(self):
        mapping = Mapping(
            self.data_dir.joinpath("martini3/naphthalene.map"), DummyOptions
        )
        self.assertEqual(1, len(mapping))
        self.assertTrue("NAPH" in mapping)
        self.assertEqual(5, len(mapping["NAPH"]))
        self.assertTrue(isinstance(mapping["NAPH"][2], VirtualMap))
        self.assertEqual(4, len(mapping["NAPH"][2].atoms))
        self.assertEqual("R1", mapping["NAPH"][2].atoms[0])
        self.assertEqual(
            1, [isinstance(bead, VirtualMap) for bead in mapping["NAPH"]].count(True)
        )

    def test_mapping_apply(self):
        mapping = Mapping(self.data_dir.joinpath("water.map"), DummyOptions)
        frame = Frame(self.data_dir.joinpath("water.gro"))
        cg_frame = mapping.apply(frame)

        self.assertEqual(frame.natoms / 3, cg_frame.natoms)

        cg_frame.save("water-cg.gro")

        self.assertTrue(
            filecmp.cmp(self.data_dir.joinpath("water-cg.gro"), "water-cg.gro")
        )
        os.remove("water-cg.gro")

    def test_mapping_duplicate_atom_weight(self):
        """Test that an atom can be duplicated in a bead specification to act as a weighting."""
        frame = Frame(self.data_dir.joinpath("water.gro"))

        mapping = Mapping(self.data_dir.joinpath("water.map"), DummyOptions)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(
            np.array([[0.716, 1.300, 1.198]]),
            cg_frame.residue(0).atom(0).coords,
            atol=0.0005,
        )

        mapping = Mapping(
            self.data_dir.joinpath("water_double_weight.map"), DummyOptions
        )
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(
            np.array([[0.720, 1.294, 1.195]]),
            cg_frame.residue(0).atom(0).coords,
            atol=0.0005,
        )

    def test_mapping_charges(self):
        mapping = Mapping(self.data_dir.joinpath("dppc.map"), DummyOptions)
        self.assertEqual(1, mapping["DPPC"][0].charge)
        self.assertEqual(-1, mapping["DPPC"][1].charge)

    def test_mapping_pbc(self):
        frame = Frame(self.data_dir.joinpath("pbcwater.gro"))

        mapping = Mapping(self.data_dir.joinpath("water.map"), DummyOptions)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(frame.atom(0).coords, cg_frame.atom(0).coords)

    def test_mapping_weights_geom(self):
        frame = Frame(self.data_dir.joinpath("two.gro"))

        mapping = Mapping(self.data_dir.joinpath("two.map"), DummyOptions)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(
            np.array([[1.5, 1.5, 1.5]]), cg_frame.residue(0).atom(0).coords
        )

    def test_virtual_mapping_weights_geom(self):
        frame = Frame(self.data_dir.joinpath("martini3/four.gro"))

        mapping = Mapping(self.data_dir.joinpath("martini3/four.map"), DummyOptions)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(
            np.array([[2.5, 2.5, 2.5]]), cg_frame.residue(0).atom(2).coords
        )

    def test_mapping_weights_mass(self):
        frame = Frame(self.data_dir.joinpath("two.gro"))
        options = DummyOptions()
        options.map_center = "mass"

        mapping = Mapping(
            self.data_dir.joinpath("two.map"),
            options,
            itp_filename=self.data_dir.joinpath("two.itp"),
        )
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(
            np.array([[2.0, 2.0, 2.0]]), cg_frame.residue(0).atom(0).coords
        )

    def test_virtual_mapping_weights_mass(self):
        frame = Frame(self.data_dir.joinpath("martini3/four.gro"))
        options = DummyOptions()
        options.virtual_map_center = "mass"

        mapping = Mapping(
            self.data_dir.joinpath("martini3/four.map"),
            options,
            itp_filename=self.data_dir.joinpath("martini3/four.itp"),
        )
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(
            np.array([[3.0, 3.0, 3.0]]), cg_frame.residue(0).atom(2).coords
        )

    def test_mapping_weights_guess_mass(self):
        frame = Frame(self.data_dir.joinpath("two.gro"))
        options = DummyOptions()
        options.map_center = "mass"

        mapping = Mapping(self.data_dir.joinpath("two.map"), options)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(
            np.array([[1.922575, 1.922575, 1.922575]], dtype=np.float32),
            cg_frame.residue(0).atom(0).coords,
            rtol=0.01,
        )

    def test_virtual_mapping_weights_guess_mass(self):
        frame = Frame(self.data_dir.joinpath("martini3/four.gro"))
        options = DummyOptions()
        options.virtual_map_center = "mass"

        mapping = Mapping(self.data_dir.joinpath("martini3/four.map"), options)
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(
            np.array([[2.83337, 2.83337, 2.83337]], dtype=np.float32),
            cg_frame.residue(0).atom(2).coords,
            rtol=0.01,
        )

    def test_mapping_weights_first(self):
        frame = Frame(self.data_dir.joinpath("two.gro"))
        options = DummyOptions()
        options.map_center = "first"

        mapping = Mapping(
            self.data_dir.joinpath("two.map"),
            options,
            itp_filename=self.data_dir.joinpath("two.itp"),
        )
        cg_frame = mapping.apply(frame)

        np.testing.assert_allclose(
            np.array([[1.0, 1.0, 1.0]]), cg_frame.residue(0).atom(0).coords
        )

    def test_mapping_itp_multi(self):
        mapping = Mapping(
            self.data_dir.joinpath("membrane/membrane.map"),
            DummyOptions,
            itp_filename=self.data_dir.joinpath("membrane/membrane.top"),
        )
        self.assertAlmostEqual(-1.2, mapping["POPE"][0].charge, delta=0.0001)
        self.assertAlmostEqual(0, mapping["POPG"][0].charge, delta=0.0001)

        self.assertAlmostEqual(94.9716, mapping["POPE"][0].mass, delta=0.0001)

    def test_empty_result_frame(self):
        """Test that no error is raised when an empty output frame is encountered.

        This is instead checked and logged as a warning in the outer script.
        """
        frame = Frame(self.data_dir.joinpath("sugar.gro"))
        mapping = Mapping(self.data_dir.joinpath("water.map"), DummyOptions)

        cg_frame = mapping.apply(frame)

        self.assertEqual(0, cg_frame.natoms)
        self.assertEqual((1, 0, 3), cg_frame._trajectory.xyz.shape)
