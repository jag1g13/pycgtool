import unittest

from pycgtool.parsers.cfg import DuplicateSectionError, NoSectionError
from pycgtool.parsers import ITP


class TestParsersITP(unittest.TestCase):
    def test_empty_itp(self):
        itp = ITP()
        self.assertEqual(0, len(itp))

    def test_itp_with(self):
        with ITP("test/data/membrane/membrane.top"):
            pass

    def test_itp_duplicate_error(self):
        with self.assertRaises(DuplicateSectionError):
            ITP("test/data/membrane/twice.itp")

    def test_itp_include(self):
        with ITP("test/data/membrane/membrane.top") as itp:
            self.assertIn('CDL2', itp)
            self.assertIn('POPE', itp)

    def test_read_section(self):
        pii_line = ["1", "CC3161", "1", "PII", "C13", "1", "0.1400", "12.0110"]

        with ITP("test/data/membrane/membrane.top") as itp:
            self.assertEqual(len(itp["PII"]["atoms"]), 283)
            self.assertListEqual(list(itp["PII"]["atoms"][0]), pii_line)

    def test_double_sections(self):
        with ITP("test/data/two_double_section.itp") as itp:
            self.assertEqual(len(itp["TWO"]["atoms"]), 2)

    def test_unknown_sections(self):
        with ITP("test/data/two_unknown_section.itp") as itp:
            self.assertEqual(len(itp["TWO"]["atoms"]), 2)

    def test_error_no_sections(self):
        with self.assertRaises(NoSectionError):
            ITP("test/data/nosections.cfg")
