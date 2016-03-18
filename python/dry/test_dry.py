import unittest
import logging
import setup_GCToo_logger as setup_logger

# import argparse
import in_out.parse_gctoo as parse_gctoo
import numpy as np
import pandas as pd
import dry

import in_out.gct as gct

logger = logging.getLogger(setup_logger.LOGGER_NAME)
PSP_config_path = "~/.PSP_config"
test_p100_path = "/Users/lev/code/PSP/python/functional_tests/p100_prm_plate29_3H.gct"
test_p100_processed_path = "/Users/lev/code/PSP/python/functional_tests/p100_prm_plate29_3H_processed.gct"

# N.B. e_out is expected_output.
class TestDry(unittest.TestCase):
    def test_main(self):
        args_string = ("{} -out_path ./ " +
                       "-out_name testing_testing " +
                       "-PSP_config_path {}").format(test_p100_path, PSP_config_path)
        args = dry.build_parser().parse_args(args_string.split())
        e_gct = parse_gctoo.parse(test_p100_processed_path)

        out_gct = dry.main(args)
        # self.assertEqual(e_gct, out_gct, "The expected and actual processed gcts are not the same.")
        pass

    def test_read_gct_and_check_provenance_code(self):
        e_prov_code = ["PRM", "L2X"]
        e_data_df_shape = (96, 96)
        (out_gct, out_prov_code) = dry.read_gct_and_check_provenance_code(PSP_config_path, test_p100_path)
        self.assertEqual(out_prov_code, e_prov_code,
                         ("The expected provenance code is {}," +
                          "not {}.".format(e_prov_code, out_prov_code)))

        self.assertEqual(out_gct.data_df.shape, e_data_df_shape,
                         ("The expected shape of the data matrix is {}," +
                          "not {}.".format(e_data_df_shape, out_gct.data_df.shape)))


    def test_extract_prov_code(self):
        col_meta_df = pd.DataFrame.from_dict({"foo": ["a", "b", "c"],
                                              "provenance_code": ["PRM+L2X",
                                                                  "PRM+L2X",
                                                                  "PRM+L2X"]},
                                             orient='index')
        e_prov_code = ["PRM", "L2X"]
        prov_code = dry.extract_prov_code(col_meta_df)
        self.assertEqual(e_prov_code, prov_code,
                        ("The expected provenance code is {}, " +
                         "not {}.").format(e_prov_code, prov_code))

    def test_check_assay_type(self):
        assay_type = "aBc"
        allowed = ["ABC", "aBc", "d"]
        assay_ok = dry.check_assay_type(assay_type, allowed)
        self.assertTrue(assay_ok)

    def test_log_transform(self):
        in_df = pd.DataFrame(np.array([[10, -3, 1.2],
                                       [0.45, 0.2, 0],
                                       [4.5, np.nan, 0.3]], dtype=float))
        e_df = pd.DataFrame(np.array([[3.322, np.nan, 0.263],
                                      [-1.152, -2.322, np.nan],
                                      [2.170, np.nan, -1.737]], dtype=float))
        out_df = dry.log_transform(in_df, log_base=2)
        self.assertTrue(np.allclose(out_df, e_df, atol=1e-3, equal_nan=True),
                        ("\nExpected output:\n{} " +
                         "\nActual output:\n{}").format(e_df, out_df))

    def test_gcp_histone_normalize(self):
        df = pd.DataFrame([[1.1, 2.0, 3.3], [4.1, 5.8, 6.0]],
                          index=["a","b"],
                          columns=["c1","c2","c3"])
        e_df = pd.DataFrame([[3.0, 3.8, 2.7]],
                          index=["b"],
                          columns=["c1","c2","c3"])
        out_df = dry.gcp_histone_normalize(df, "a")
        self.assertTrue(out_df.shape == e_df.shape,
                        ("expected_out.shape: {} not the same " +
                         "as actual_out.shape: {}").format(e_df.shape, out_df.shape))
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
        medians = np.array([2, 0.5, 0.1, 1.1, 1.5])
        e_out = 45.62
        out = dry.distance_function(offset, values, medians)
        self.assertTrue(np.isclose(out, e_out),
                        ("\nExpected output: {} " +
                         "\nActual output: {}").format(e_out, out))

        offset = 0
        values = np.array([1, 2.5, np.nan, 4, 5])
        medians = np.array([2, 0.5, 0.1, 1.1, 1.5])
        e_out = 25.66
        out = dry.distance_function(offset, values, medians)
        self.assertTrue(np.isclose(out, e_out),
                        ("\nExpected output: {} " +
                         "\nActual output: {}").format(e_out, out))

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

    def test_calculate_distances_only(self):
        df = pd.DataFrame(np.array([[10, 3, 1.2],
                                    [0.45, 0.2, np.nan],
                                    [4.5, 4, 0.3]], dtype=float))
        e_dists = [49.27, 0.02, 16.93]
        out_dists = dry.calculate_distances_only(df)
        self.assertTrue(np.allclose(e_dists, out_dists, atol=1e-2),
                        ("The expected distances are {}, " +
                         "not {}.").format(e_dists, out_dists))

    def test_calculate_distances_and_optimize_if_needed(self):
        df = pd.DataFrame(np.array([[10, 3, 1.2],
                                    [0.45, 0.2, np.nan],
                                    [4.5, 4, 0.3]], dtype=float))

        # optim_bool = False, assay_type = GCP
        optim_bool = False
        prov_code = ["GR1", "L2X"]
        e_prov_code = prov_code
        (_, _, _, out_prov_code) = dry.calculate_distances_and_optimize_if_needed(
            df, optim_bool, (-7, 7), prov_code)
        self.assertTrue(np.array_equal(e_prov_code, out_prov_code), (
            "The expected output provenance code is {}, not {}.".format(e_prov_code, out_prov_code)))

        # optim_bool = False, assay_type = P100
        optim_bool = False
        prov_code = ["PR1", "L2X"]
        e_prov_code = prov_code
        (_, _, _, out_prov_code) = dry.calculate_distances_and_optimize_if_needed(
            df, optim_bool, (-7, 7), prov_code)
        self.assertTrue(np.array_equal(e_prov_code, out_prov_code), (
            "The expected output provenance code is {}, not {}.".format(e_prov_code, out_prov_code)))

        # optim_bool = True, assay_type = GCP
        optim_bool = True
        prov_code = ["GR1", "L2X"]
        e_prov_code = prov_code
        (_, _, _, out_prov_code) = dry.calculate_distances_and_optimize_if_needed(
            df, optim_bool, (-7, 7), prov_code)
        self.assertTrue(np.array_equal(e_prov_code, out_prov_code), (
            "The expected output provenance code is {}, not {}.".format(e_prov_code, out_prov_code)))

        # optim_bool = True, assay_type = P100
        optim_bool = True
        prov_code = ["DIA1", "L2X"]
        e_prov_code = np.append(prov_code, "LLB")
        (_, _, _, out_prov_code) = dry.calculate_distances_and_optimize_if_needed(
            df, optim_bool, (-7, 7), prov_code)
        self.assertTrue(np.array_equal(e_prov_code, out_prov_code), (
            "The expected output provenance code is {}, not {}.".format(e_prov_code, out_prov_code)))



    def test_remove_sample_outliers(self):
        df = pd.DataFrame(np.array([[10, -3, 1.2, 0.6],
                                    [0.45, 0.2, 0, 0.2],
                                    [4.5, np.nan, 0.3, 0.4]], dtype=float))
        dists = np.array([0.2, 5, 0.5, 0.4], dtype=float)
        bools = np.array([True, True, True, False], dtype=bool)
        out = dry.remove_sample_outliers(df, dists, bools,
                                         dist_sd_cutoff=1)
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
