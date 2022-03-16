import pathlib
import unittest

from pycgtool.parsers.cfg import DuplicateSectionError, NoSectionError
from pycgtool.parsers import CFG


class TestParsersCFG(unittest.TestCase):
    base_dir = pathlib.Path(__file__).absolute().parent
    data_dir = base_dir.joinpath("data")

    water_name = "SOL"
    watline = ("W", "P4", "OW", "HW1", "HW2")

    def test_cfg_with(self):
        with CFG(self.data_dir.joinpath("water.map")):
            pass

    def test_cfg_iterate_sections(self):
        with CFG(self.data_dir.joinpath("water.map")) as cfg:
            for name in cfg:
                self.assertEqual(self.water_name, name)

    def test_cfg_iterate_lines(self):
        with CFG(self.data_dir.joinpath("water.map")) as cfg:
            for name, section in cfg.items():
                self.assertEqual(self.water_name, name)
                for line in section:
                    self.assertEqual(self.watline, line)

    def test_cfg_get_section(self):
        with CFG(self.data_dir.joinpath("water.map")) as cfg:
            self.assertTrue(self.water_name in cfg)
            for line in cfg[self.water_name]:
                self.assertEqual(self.watline, line)

    def test_cfg_duplicate_error(self):
        with self.assertRaises(DuplicateSectionError):
            CFG(self.data_dir.joinpath("twice.cfg"))

    def test_include_file(self):
        with CFG(self.data_dir.joinpath("martini.map")) as cfg:
            self.assertTrue("DOPC" in cfg)
            self.assertTrue("GLY" in cfg)

    def test_error_no_sections(self):
        with self.assertRaises(NoSectionError):
            CFG(self.data_dir.joinpath("nosections.cfg"))


if __name__ == "__main__":
    unittest.main()
