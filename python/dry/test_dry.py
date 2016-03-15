import unittest
import logging
import setup_GCToo_logger as setup_logger
import numpy as np
import pandas as pd
import dry
import os
import in_out.gct as gct

logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestDry(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.func_tests_path = "/Users/lev/code/PSP/python/functional_tests"

    def test_extract_prov_code(self):
        col_metadata_df = pd.DataFrame(np.array([["a", "b", "c"], ["PRM+L2X", "PRM+L2X", "PRM+L2X"]]),
                                       index=["foo", "provenance_code"])
        expected_prov_code = ["PRM", "L2X"]
        prov_code = dry.extract_prov_code(col_metadata_df)
        self.assertTrue(expected_prov_code == prov_code,
                        ("The expected provenance code is {}, " +
                         "but actual provenance code is {}.").format(expected_prov_code, prov_code))

    def test_update_prov_code(self):
        existing_prov_code = np.array(["PR1"])
        new_entry = "L2X"
        expected_out = np.array(["PR1", "L2X"])
        actual_out = dry.update_prov_code(new_entry, existing_prov_code)
        self.assertTrue(np.array_equal(actual_out, expected_out),
                        "Expected output: {}, actual output: {}".format(expected_out, actual_out))

    def test_log_transform(self):
        data_df = pd.DataFrame(np.array([[10, -3, 1.2],
                                         [0.45, 0.2, 0],
                                         [4.5, np.nan, 0.3]], dtype=float))
        expected_out = pd.DataFrame(np.array([[3.322, np.nan, 0.263],
                                              [-1.152, -2.322, np.nan],
                                              [2.170, np.nan, -1.737]], dtype=float))
        actual_out = dry.log_transform(data_df)
        self.assertTrue(np.allclose(actual_out, expected_out, atol=1e-3, equal_nan=True),
                        "\nExpected output:\n{} \nActual output:\n{}".format(expected_out, actual_out))

    def test_manual_probe_rejection(self):
        row_metadata_df = pd.DataFrame({"foo": ["a", "b", "c"],
                                        "pr_probe_suitability_manual": ["TRUE", "TRUE", "FALSE"],
                                        "bar": ["d", "e", "f"]})
        data_df = pd.DataFrame(np.array([[10, -3, 1.2],
                                         [0.45, 0.2, 0],
                                         [4.5, np.nan, 0.3]], dtype=float))
        expected_meta_out = pd.DataFrame({"foo": ["a", "b"],
                                          "pr_probe_suitability_manual": ["TRUE", "TRUE"],
                                          "bar": ["d", "e"]})
        expected_data_out = pd.DataFrame(np.array([[10, -3, 1.2], [0.45, 0.2, 0]], dtype=float))
        (actual_data_out, actual_meta_out) = dry.manual_probe_rejection(data_df, row_metadata_df)

        self.assertTrue(np.array_equal(actual_meta_out, expected_meta_out),
                        ("\nExpected meta_out:\n{} " +
                         "\nActual meta_out:\n{}").format(expected_meta_out, actual_meta_out))
        self.assertTrue(np.allclose(actual_data_out, expected_data_out, atol=1e-3, equal_nan=True),
                        ("\nExpected data_out:\n{} " +
                         "\nActual data_out:\n{}").format(expected_data_out, actual_data_out))

    def test_filter_samples(self):
        df = pd.DataFrame(np.array([[0.5, 0.2, 0.1, 0.25],
                                    [np.nan, 0.45, 0.2, -0.1],
                                    [np.nan, 0.02, np.nan, 0.3]], dtype=float))
        expected_out = pd.DataFrame(np.array([[0.2, 0.1, 0.25],
                                              [0.45, 0.2, -0.1],
                                              [0.02, np.nan, 0.3]], dtype=float))
        actual_out = dry.filter_samples(df, sample_pct_cutoff=0.4)
        self.assertTrue(actual_out.shape == expected_out.shape,
                        ("expected_out.shape: {} not the same " +
                         "as actual_out.shape: {}").format(expected_out.shape, actual_out.shape))

        self.assertTrue(np.allclose(actual_out, expected_out, equal_nan=True),
                        "\nExpected output:\n{} \nActual output:\n{}".format(expected_out, actual_out))

    def test_filter_probes(self):
        df = pd.DataFrame(np.array([[10, 0.2, 0.1, 0.25],
                                    [np.nan, 0.45, 0.2, -0.1],
                                    [np.nan, 0.02, np.nan, 0.3]], dtype=float))
        expected_out = pd.DataFrame(np.array([[np.nan, 0.45, 0.2, -0.1]], dtype=float))
        actual_out = dry.filter_probes(df, probe_pct_cutoff=0.3, probe_sd_cutoff=3)
        self.assertTrue(actual_out.shape == expected_out.shape,
                        ("expected_out.shape: {} not the same " +
                         "as actual_out.shape: {}").format(expected_out.shape, actual_out.shape))
        self.assertTrue(np.allclose(actual_out, expected_out, equal_nan=True),
                        "\nExpected output:\n{} \nActual output:\n{}".format(expected_out, actual_out))

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
                                                          [0.08, -1.16, 0.40]],
                                                         dtype=float))
        expected_distances = np.array([36.62, 12.04, 0.06], dtype=float)
        (actual_df_with_offsets,
         actual_distances,
         success_bools) = dry.optimize_sample_balance(df)
        self.assertTrue(np.allclose(actual_distances, expected_distances, atol=1e-2),
                        ("\nExpected output:\n{} " +
                         "\nActual output:\n{}").format(expected_distances,
                                                        actual_distances))
        self.assertTrue(np.allclose(actual_df_with_offsets,
                                    expected_df_with_offsets, atol=1e-2),
                        ("\nExpected output:\n{} " +
                         "\nActual output:\n{}").format(expected_df_with_offsets,
                                                        actual_df_with_offsets))
        self.assertTrue(~all(success_bools),
                        "All samples should have converged.")

    def test_remove_sample_outliers(self):
        df = pd.DataFrame(np.array([[10, -3, 1.2, 0.6],
                                    [0.45, 0.2, 0, 0.2],
                                    [4.5, np.nan, 0.3, 0.4]], dtype=float))
        dists = np.array([0.2, 5, 0.5, 0.4], dtype=float)
        bools = np.array([True, True, True, False], dtype=bool)
        actual_out = dry.remove_sample_outliers(df, dists, bools,
                                                sd_sample_outlier_cutoff=1)
        e_out = df.iloc[:, [0, 2]]
        self.assertTrue(actual_out.shape == e_out.shape,
                        ("expected_out.shape: {} not the same " +
                         "as actual_out.shape: {}").format(e_out.shape,
                                                           actual_out.shape))



    def test_main(self):
        dry.main()


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()