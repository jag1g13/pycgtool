import unittest
import subprocess
import os
import logging

import numpy as np

from simpletraj.trajectory import XtcTrajectory

try:
    import mdtraj
    mdtraj_present = True
except ImportError:
    mdtraj_present = False

from pycgtool.interface import Options
from pycgtool.util import cmp_whitespace_float
from pycgtool.pycgtool import main, map_only


class Args:
    def __init__(self, name, map=True, bnd=True):
        self.gro = os.path.join("test/data", name+".gro")
        self.xtc = os.path.join("test/data", name+".xtc")
        self.map = os.path.join("test/data", name+".map") if map else None
        self.bnd = os.path.join("test/data", name+".bnd") if bnd else None
        self.begin = 0
        self.end = -1
        self.quiet = True


class PycgtoolTest(unittest.TestCase):
    config = Options([("output_name", "out"),
                      ("output", "gro"),
                      ("output_xtc", True),
                      ("map_only", False),
                      ("map_center", "geom"),
                      ("constr_threshold", 100000),
                      ("dump_measurements", False),
                      ("dump_n_values", 10000),
                      ("output_forcefield", False),
                      ("temperature", 310),
                      ("default_fc", False),
                      ("generate_angles", True),
                      ("generate_dihedrals", False)])

    def test_run_help(self):
        path = os.path.dirname(os.path.dirname(__file__))
        self.assertEqual(0, subprocess.check_call([os.path.join(path, "pycgtool.py"), "-h"], stdout=subprocess.PIPE))

    @unittest.skipIf(not mdtraj_present, "MDTRAJ or Scipy not present")
    def test_map_only(self):
        logging.disable(logging.WARNING)
        map_only(Args("sugar"), self.config)
        logging.disable(logging.NOTSET)

        xtc = XtcTrajectory("out.xtc")
        xtc_ref = XtcTrajectory("test/data/sugar_out.xtc")
        self.assertEqual(xtc_ref.numframes, xtc.numframes)

        for i in range(xtc_ref.numframes):
            xtc.get_frame(i)
            xtc_ref.get_frame(i)
            np.testing.assert_array_almost_equal(xtc_ref.box, xtc.box, decimal=3)
            np.testing.assert_array_almost_equal(xtc_ref.x, xtc.x, decimal=3)

    def test_full(self):
        path = os.path.dirname(os.path.dirname(__file__))
        self.assertEqual(0, subprocess.check_call([os.path.join(path, "pycgtool.py"),
                                                   "-g", "test/data/sugar.gro",
                                                   "-x", "test/data/sugar.xtc",
                                                   "-m", "test/data/sugar_only.map",
                                                   "-b", "test/data/sugar.bnd",
                                                   ], stdout=subprocess.PIPE, stderr=subprocess.PIPE))
        self.assertTrue(cmp_whitespace_float("out.itp", "test/data/sugar_out.itp", float_rel_error=0.001))
        self.assertTrue(cmp_whitespace_float("out.gro", "test/data/sugar_out.gro", float_rel_error=0.001))
    # TODO more tests


if __name__ == '__main__':
    unittest.main()
