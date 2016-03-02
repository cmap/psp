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
        actual_out = dry.filter_samples(df, sample_pct_cutoff=0.4)
        self.assertTrue(np.allclose(actual_out, expected_out, equal_nan=True),
                        "\nExpected output:\n{} \nActual output:\n{}".format(expected_out, actual_out))

    def test_filter_probes(self):
        df = pd.DataFrame(np.array([[10, 0.2, 0.1, 0.25],
                                    [np.nan, 0.45, 0.2, -0.1],
                                    [np.nan, 0.02, np.nan, 0.3]], dtype=float))
        expected_out = np.array([np.nan, 0.45, 0.2, -0.1], dtype=float)
        actual_out = dry.filter_probes(df, probe_pct_cutoff=0.3, probe_sd_cutoff=3)
        self.assertTrue(np.allclose(actual_out.values, expected_out, equal_nan=True),
                        "\nExpected output:\n{} \nActual output:\n{}".format(expected_out, actual_out.values))

    def test_function_to_optimize(self):
        offset = 0.5
        values = np.array([1, 2.5, 3, 4, 5])
        means = np.array([2, 0.5, 0.1, 1.1, 1.5])
        expected_out = 45.62
        actual_out = dry.function_to_optimize(offset, values, means)
        self.assertTrue(np.isclose(actual_out, expected_out),
                        "\nExpected output: {} \nActual output: {}".format(expected_out, actual_out))

    def test_optimize_sample_balance(self):
        df = pd.DataFrame(np.array([[10, -3, 1.2],
                                    [0.45, 0.2, -0.1],
                                    [4.5, -4, 0.3]], dtype=float))
        expected_df_with_offsets = pd.DataFrame(np.array([[5.58, -0.16, 1.30],
                                                          [-3.96, 3.03, 0],
                                                          [0.08, -1.16, 0.40]], dtype=float))
        expected_distances = np.array([36.62,  12.04,  0.06], dtype=float)
        (actual_df_with_offsets, actual_distances) = dry.optimize_sample_balance(df)
        self.assertTrue(np.allclose(actual_distances, expected_distances, atol=1e-2),
                        "\nExpected output:\n{} \nActual output:\n{}".format(expected_distances,
                                                                             actual_distances))
        self.assertTrue(np.allclose(actual_df_with_offsets, expected_df_with_offsets, atol=1e-2),
                        "\nExpected output:\n{} \nActual output:\n{}".format(expected_df_with_offsets,
                                                                             actual_df_with_offsets))

    def test_remove_sample_outliers(self):
        # Haven't been able to create an example of an outlier sample
        df = pd.DataFrame(np.array([[2.4, -0.1, 5.41],
                                    [-1, 0.3, 0.20],
                                    [1.5, -1.2, 0.40]], dtype=float))
        expected_df = df
        (_, distances) = dry.optimize_sample_balance(df)
        dry.remove_sample_outliers(df, distances, sd_sample_outlier_cutoff=2)

    def test_main(self):
        dry.main()


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()