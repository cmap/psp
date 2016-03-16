import unittest
import logging
import setup_GCToo_logger as setup_logger
import numpy as np
import pandas as pd
import dry
import os
import in_out.gct as gct

logger = logging.getLogger(setup_logger.LOGGER_NAME)


### TO-DO(lev): Test all the do_* methods.

# N.B. e_out is expected_output.
class TestDry(unittest.TestCase):
    def test_main(self):
        # input = unprocessed
        # check that output = processed gct
        pass

    def test_extract_prov_code(self):
        col_meta_df = pd.DataFrame.from_dict({"foo": ["a", "b", "c"],
                                              "provenance_code": ["PRM+L2X",
                                                                  "PRM+L2X",
                                                                  "PRM+L2X"]},
                                             orient='index')
        e_prov_code = ["PRM", "L2X"]
        prov_code = dry.extract_prov_code(col_meta_df)
        self.assertTrue(e_prov_code == prov_code,
                        ("The expected provenance code is {}, " +
                         "not {}.").format(e_prov_code, prov_code))

    def test_check_assay_type(self):
        assay_type = "aBc"
        allowed = ["ABC", "aBc", "d"]
        assay_ok = dry.check_assay_type(assay_type, allowed)
        self.assertTrue(assay_ok)

    def test_do_log_transform_if_needed(self):
        pass

    def test_log_transform(self):
        in_df = pd.DataFrame(np.array([[10, -3, 1.2],
                                       [0.45, 0.2, 0],
                                       [4.5, np.nan, 0.3]], dtype=float))
        e_df = pd.DataFrame(np.array([[3.322, np.nan, 0.263],
                                      [-1.152, -2.322, np.nan],
                                      [2.170, np.nan, -1.737]], dtype=float))
        out_df = dry.log_transform(in_df)
        self.assertTrue(np.allclose(out_df, e_df, atol=1e-3, equal_nan=True),
                        ("\nExpected output:\n{} " +
                         "\nActual output:\n{}").format(e_df, out_df))

    def test_filter_samples_by_nan(self):
        df = pd.DataFrame(np.array([[0.5, 0.2, 0.1, 0.25],
                                    [np.nan, 0.45, 0.2, -0.1],
                                    [np.nan, 0.02, np.nan, 0.3]], dtype=float))
        e_out = pd.DataFrame(np.array([[0.2, 0.1, 0.25],
                                       [0.45, 0.2, -0.1],
                                       [0.02, np.nan, 0.3]], dtype=float))
        out = dry.filter_samples_by_nan(df, sample_nan_thresh=0.6)
        self.assertTrue(out.shape == e_out.shape,
                        ("expected_out.shape: {} not the same " +
                         "as actual_out.shape: {}").format(e_out.shape, out.shape))

        self.assertTrue(np.allclose(out, e_out, equal_nan=True),
                        ("\nExpected output:\n{} " +
                         "\nActual output:\n{}").format(e_out, out))

    def test_manual_probe_rejection(self):
        row_meta_df = pd.DataFrame({"foo": ["a", "b", "c"],
                                    "pr_probe_suitability_manual": ["TRUE", "TRUE", "FALSE"],
                                    "bar": ["d", "e", "f"]})
        data_df = pd.DataFrame(np.array([[10, -3, 1.2],
                                         [0.45, 0.2, 0],
                                         [4.5, np.nan, 0.3]], dtype=float))
        e_df = pd.DataFrame(np.array([[10, -3, 1.2], [0.45, 0.2, 0]], dtype=float))
        out_df = dry.manual_probe_rejection(data_df, row_meta_df)

        self.assertTrue(np.allclose(out_df, e_df, atol=1e-3, equal_nan=True),
                        ("\nExpected df:\n{} " +
                         "\nActual df:\n{}").format(e_df, out_df))

    def test_filter_probes_by_nan_and_sd(self):
        df = pd.DataFrame(np.array([[10, 0.2, 0.1, 0.25],
                                    [np.nan, 0.45, 0.2, -0.1],
                                    [np.nan, 0.02, np.nan, 0.3]], dtype=float))
        e_out = df.iloc[[1], :]
        out = dry.filter_probes_by_nan_and_sd(df, probe_nan_thresh=0.6, probe_sd_cutoff=3)
        self.assertTrue(out.shape == e_out.shape,
                        ("expected_out.shape: {} not the same " +
                         "as actual_out.shape: {}").format(e_out.shape, out.shape))
        self.assertTrue(np.allclose(out, e_out, equal_nan=True),
                        ("\nExpected output:\n{} " +
                         "\nActual output:\n{}").format(e_out, out))

    def test_distance_function(self):
        offset = 0.5
        values = np.array([1, 2.5, 3, 4, 5])
        means = np.array([2, 0.5, 0.1, 1.1, 1.5])
        e_out = 45.62
        out = dry.distance_function(offset, values, means)
        self.assertTrue(np.isclose(out, e_out),
                        ("\nExpected output: {} " +
                         "\nActual output: {}").format(e_out, out))

    def test_calculate_distances_and_optimize_if_needed(self):
        pass

    def test_calculate_distances_and_optimize(self):
        df = pd.DataFrame(np.array([[10, -3, 1.2],
                                    [0.45, 0.2, -0.1],
                                    [4.5, -4, 0.3]], dtype=float))
        e_df = pd.DataFrame(np.array([[5.58, -0.16, 1.30],
                                      [-3.96, 3.03, 0],
                                      [0.08, -1.16, 0.40]],
                                     dtype=float))
        e_dists = np.array([36.62, 12.04, 0.06], dtype=float)
        (out_df, dists, success_bools) = dry.calculate_distances_and_optimize(df, (-7,7))
        self.assertTrue(np.allclose(dists, e_dists, atol=1e-2),
                        ("\nExpected distances:\n{} " +
                         "\nActual distances:\n{}").format(e_dists, dists))
        self.assertTrue(np.allclose(out_df, e_df, atol=1e-2),
                        ("\nExpected out_df:\n{} " +
                         "\nActual out_df:\n{}").format(e_df, out_df))
        self.assertTrue(~all(success_bools),
                        "All samples should have converged.")

    def test_calculate_distances(self):
        pass

    def test_distance_function(self):
        pass

    def test_remove_sample_outliers(self):
        df = pd.DataFrame(np.array([[10, -3, 1.2, 0.6],
                                    [0.45, 0.2, 0, 0.2],
                                    [4.5, np.nan, 0.3, 0.4]], dtype=float))
        dists = np.array([0.2, 5, 0.5, 0.4], dtype=float)
        bools = np.array([True, True, True, False], dtype=bool)
        out = dry.remove_sample_outliers(df, dists, bools,
                                         sd_sample_outlier_cutoff=1)
        e_out = df.iloc[:, [0, 2]]
        self.assertTrue(out.shape == e_out.shape, (
            "expected_out.shape: {} not the same " +
            "as actual_out.shape: {}").format(e_out.shape, out.shape))


    def test_row_median_normalize(self):
        df = pd.DataFrame(np.array([[10, -3, 1.2, 0.6],
                                    [0.45, 0.2, 0, 0.2],
                                    [4.5, np.nan, 0.3, 0.4]], dtype=float))
        e_df = pd.DataFrame(np.array([[9.1, -3.9, 0.3, -0.3],
                                    [0.25, 0, -0.2, 0],
                                    [4.1, np.nan, -0.1, 0]], dtype=float))
        out_df = dry.row_median_normalize(df)
        self.assertTrue(np.allclose(out_df, e_df, atol=1e-2, equal_nan=True),
                        ("\nExpected out_df:\n{} " +
                         "\nActual out_df:\n{}").format(e_df, out_df))


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
