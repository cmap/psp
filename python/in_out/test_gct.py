import unittest
import logging
import utils.setup_logger as setup_logger
import ConfigParser

import os
import gct
import numpy as np

logger = logging.getLogger(setup_logger.LOGGER_NAME)
configParser = ConfigParser.RawConfigParser()

config_path = os.path.expanduser("~/.PSP_config")
configParser.read(config_path)
functional_tests_path = configParser.get("tests", "functional_tests_path")


class TestGct(unittest.TestCase):
    def test_strlist2floatlist(self):
        # Case A: no NaNs
        list_wo_nan = ["5.674", "-7.263", "0.152"]
        expected_out = [5.674, -7.263, 0.152]
        actual_out = gct.GCT.strlist2floatlist(list_wo_nan)

        # Convert outputs to ndarrays in order to use np.allclose
        expected_out_np = np.array(expected_out, dtype=np.float)
        actual_out_np = np.array(actual_out, dtype=np.float)
        self.assertTrue(np.allclose(expected_out_np, actual_out_np, equal_nan=True),
                        "Expected output: {}, Actual output: {}".format(expected_out, actual_out))

        # Case B: NaNs present
        list_with_nan = ["NaN", "5.457", "NA", "NULL", "#N/A"]
        expected_out = [None, 5.457, None, None, None]
        actual_out = gct.GCT.strlist2floatlist(list_with_nan)

        # Convert outputs to ndarrays in order to use np.allclose
        expected_out_np = np.array(expected_out, dtype=np.float)
        actual_out_np = np.array(actual_out, dtype=np.float)
        self.assertTrue(np.allclose(expected_out_np, actual_out_np, equal_nan=True),
                        "Expected output: {}, Actual output: {}".format(expected_out, actual_out))

    def test_read_gct(self):
        # gct.py cannot read a subset of gct files, but it can do this for gctx
        gct_path = os.path.join(functional_tests_path, "gct_v13.gct")
        logger.debug("gct_path: {}".format(gct_path))
        gct_obj = gct.GCT(gct_path)
        gct_obj.read()

    def test_read_p100(self):
        gct_path = os.path.join(functional_tests_path, "test_p100.gct")
        logger.debug("gct_path: {}".format(gct_path))
        gct_obj = gct.GCT(gct_path)
        gct_obj.read()

    def test_read_GCP(self):
        gct_path = os.path.join(functional_tests_path, "test_GCP.gct")
        logger.debug("gct_path: {}".format(gct_path))
        gct_obj = gct.GCT(gct_path)
        gct_obj.read()

    def test_read_gctx(self):
        # Very possibly will fail if gctx contains NaNs.
        gct_path = os.path.join(functional_tests_path, "gct.gctx")
        logger.debug("gct_path: {}".format(gct_path))
        gct_obj = gct.GCT(gct_path)
        gct_obj.read(row_inds=range(3), col_inds=range(4))


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()