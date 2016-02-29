import unittest
import logging
import utils.setup_logger as setup_logger
import numpy as np
import pandas as pd
import dry
import os
import in_out.gct as gct

logger = logging.getLogger(setup_logger.LOGGER_NAME)
functional_tests_path = "/Users/lev/code/PSP/python/functional_tests"

class TestDry(unittest.TestCase):
    def test_determine_assay_type(self):
        prov_code = "GR1"
        expected_out = np.array("GR1", dtype=np.str)
        actual_out = dry.determine_assay_type(prov_code)
        self.assertTrue(np.array_equal(actual_out, expected_out),
                        "Expected output: {}, Actual output: {}".format(expected_out, actual_out))

    def test_update_prov_code(self):
        existing_prov_code = np.array(["PR1"])
        new_entry = "L2X"
        expected_out = np.array(["PR1", "L2X"])
        actual_out = dry.update_prov_code(new_entry, existing_prov_code)
        self.assertTrue(np.array_equal(actual_out, expected_out),
                        "Expected output: {}, Actual output: {}".format(expected_out, actual_out))

    def test_filter_samples(self):
        df = pd.DataFrame(np.array([[0.5, 0.2, 0.1, 0.25],
                                    [np.nan, 0.45, 0.2, -0.1],
                                    [np.nan, 0.02, np.nan, 0.3]], dtype=float))
        expected_out = np.array([[0.2, 0.1, 0.25],
                                 [0.45, 0.2, -0.1],
                                 [0.02, np.nan, 0.3]], dtype=float)
        actual_out = dry.filter_samples(df, sample_pct_cutoff = 0.3)
        self.assertTrue(np.allclose(actual_out, expected_out, equal_nan=True),
                        "\nExpected output: {}, \nActual output: {}".format(expected_out, actual_out))

    def test_filter_probes(self):
        # filter_probes(df, probe_pct_cutoff, probe_sd_cutoff)
        # Identify rows manually labeled for rejection
        # Return df (potentially of different size than original df)
        pass

    def test_main(self):
        dry.main()


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()