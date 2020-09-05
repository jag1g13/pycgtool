import unittest

from pycgtool.parsers.cfg import DuplicateSectionError, NoSectionError
from pycgtool.parsers import ITP


class TestParsersITP(unittest.TestCase):

    def test_itp_with(self):
        with ITP("test/data/membrane/membrane.top"):
            pass

    def test_itp_duplicate_error(self):
        with self.assertRaises(DuplicateSectionError):
            ITP("test/data/membrane/twice.itp")

    def test_itp_include(self):
        with ITP("test/data/membrane/membrane.top") as itp:
            itp["CDL2"]
            itp["POPE"]

    def test_read_section(self):
        pii_line = ["1", "CC3161", "1", "PII", "C13", "1", "0.1400", "12.0110"]
        with ITP("test/data/membrane/membrane.top") as itp:
            self.assertEqual(len(itp["PII"]["atoms"]), 283)
            self.assertListEqual(list(itp["PII"]["atoms"][0]), pii_line)