import unittest

from pycgtool.parsers.cfg import DuplicateSectionError, NoSectionError
from pycgtool.parsers import CFG


class TestParsersCFG(unittest.TestCase):
    watline = ("W", "P4", "OW", "HW1", "HW2")

    def test_cfg_with(self):
        with CFG("test/data/water.map"):
            pass

    def test_cfg_iterate_sections(self):
        with CFG("test/data/water.map") as cfg:
            for name in cfg:
                self.assertEqual("SOL", name)

    def test_cfg_iterate_lines(self):
        with CFG("test/data/water.map") as cfg:
            for name, section in cfg.items():
                self.assertEqual("SOL", name)
                for line in section:
                    self.assertEqual(self.watline, line)

    def test_cfg_get_section(self):
        with CFG("test/data/water.map") as cfg:
            self.assertTrue("SOL" in cfg)
            for line in cfg["SOL"]:
                self.assertEqual(self.watline, line)

    def test_cfg_duplicate_error(self):
        with self.assertRaises(DuplicateSectionError):
            CFG("test/data/twice.cfg")

    def test_include_file(self):
        with CFG("test/data/martini.map") as cfg:
            self.assertTrue("DOPC" in cfg)
            self.assertTrue("GLY" in cfg)

    def test_error_no_sections(self):
        with self.assertRaises(NoSectionError):
            CFG("test/data/nosections.cfg")


if __name__ == '__main__':
    unittest.main()
