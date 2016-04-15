import unittest
import logging
import setup_GCToo_logger as setup_logger
import ConfigParser
import os
import sys

import in_out.parse_gctoo as parse_gctoo
import numpy as np
import pandas as pd
import dry

import in_out.gct as gct

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Set default config filepath
DEFAULT_PSP_CONFIG_PATH = "~/.PSP_config"

# N.B. e_out is expected_output.
class TestDry(unittest.TestCase):

    def test_p100_main(self):
        # Get input asset from config file
        configParser = ConfigParser.RawConfigParser()
        configParser.read(os.path.expanduser(DEFAULT_PSP_CONFIG_PATH))
        input_gct_path = configParser.get("tests", "input_p100_gct_path")

        # Write output to the functional tests directory
        functional_tests_dir = configParser.get("tests", "functional_tests_dir")

        # Specify output file name
        OUT_NAME = "test_dry_p100_output.gct"

        args_string = ("{} -out_path {} " +
                       "-out_name {} " +
                       "-PSP_config_path {} " +
                       "-sample_nan_thresh {} " +
                       "-probe_nan_thresh {} " +
                       "-probe_sd_cutoff {} " +
                       "-dist_sd_cutoff {} " +
                       "-optim -v").format(input_gct_path, functional_tests_dir,
                                           OUT_NAME, DEFAULT_PSP_CONFIG_PATH,
                                           0.8, 0.9, 3, 3)
        args = dry.build_parser().parse_args(args_string.split())
        out_gct = dry.main(args)

        # Read in output gct created by JJ's code
        e_output_gct_path = configParser.get("tests", "output_p100_gct_path")
        e_gct = parse_gctoo.parse(e_output_gct_path)

        self.assertTrue(np.allclose(out_gct.data_df, e_gct.data_df, atol=1e-1, equal_nan=True),
                        "The data_df that was returned is wrong.")
        self.assertTrue(np.array_equal(e_gct.row_metadata_df, out_gct.row_metadata_df),
                         "The row_metadata_df that was returned is wrong.")
        # Column metadata should have one more header
        self.assertEqual(e_gct.col_metadata_df.shape[1] + 1, out_gct.col_metadata_df.shape[1],
                         ("Actual col_metadata_df should have one more header" +
                          "than e_col_metadata_df.\n" +
                          "e_gct.col_metadata_df.shape: {}, " +
                          "out_gct.col_metadata_df.shape: {}").format(e_gct.col_metadata_df.shape,
                                                                      out_gct.col_metadata_df.shape[1]))
        # Check that column metadata is correct for all but the new header
        self.assertTrue(np.array_equal(e_gct.col_metadata_df, out_gct.col_metadata_df.iloc[:, :-1]),
                        "The col_metadata_df that was returned is wrong.")

        logger.info("Mean difference between gcts: {}".format(np.nanmean(e_gct.data_df.values - out_gct.data_df.values)))
        logger.info("Mean absolute difference between gcts: {}".format(np.nanmean(np.absolute(e_gct.data_df.values - out_gct.data_df.values))))

        # Clean up
        os.remove(os.path.join(functional_tests_dir, OUT_NAME))

    def test_GCP_main(self):
        # Get input asset from config file
        configParser = ConfigParser.RawConfigParser()
        configParser.read(os.path.expanduser(DEFAULT_PSP_CONFIG_PATH))
        input_gct_path = configParser.get("tests", "input_GCP_gct_path")

        # Write output to the functional tests directory
        functional_tests_dir = configParser.get("tests", "functional_tests_dir")

        # Specify output file name
        OUT_NAME = "test_dry_GCP_output.gct"

        args_string = ("{} -out_path {} " +
                       "-out_name {} " +
                       "-PSP_config_path {} " +
                       "-sample_nan_thresh {} " +
                       "-probe_nan_thresh {} " +
                       "-probe_sd_cutoff {} " +
                       "-optim -v").format(input_gct_path, functional_tests_dir,
                                           OUT_NAME, DEFAULT_PSP_CONFIG_PATH,
                                           0.5, 0.5, 4)
        args = dry.build_parser().parse_args(args_string.split())
        out_gct = dry.main(args)

        # Read in output gct created by JJ's code
        e_output_gct_path = configParser.get("tests", "output_GCP_gct_path")
        e_gct = parse_gctoo.parse(e_output_gct_path)

        self.assertTrue(np.allclose(out_gct.data_df, e_gct.data_df, atol=1e-1, equal_nan=True),
                        "The data_df that was returned is wrong.")
        self.assertTrue(np.array_equal(e_gct.row_metadata_df, out_gct.row_metadata_df),
                        "The row_metadata_df that was returned is wrong.")
        self.assertTrue(np.array_equal(e_gct.col_metadata_df, out_gct.col_metadata_df),
                        "The col_metadata_df that was returned is wrong.")

        # self.assertEqual(e_gct.col_metadata_df.shape[1], out_gct.col_metadata_df.shape[1],
        #                  ("Actual col_metadata_df should have the same number of headers " +
        #                   "as e_col_metadata_df.\n" +
        #                   "e_gct.col_metadata_df.shape: {}, " +
        #                   "out_gct.col_metadata_df.shape: {}").format(e_gct.col_metadata_df.shape,
        #                                                               out_gct.col_metadata_df.shape))

        logger.info("Mean difference between gcts: {}".format(np.nanmean(e_gct.data_df.values - out_gct.data_df.values)))
        logger.info("Mean absolute difference between gcts: {}".format(np.nanmean(np.absolute(e_gct.data_df.values - out_gct.data_df.values))))

        # Clean up
        os.remove(os.path.join(functional_tests_dir, OUT_NAME))

    def test_p100_subset_main(self):
        # Get input asset from config file
        configParser = ConfigParser.RawConfigParser()
        configParser.read(os.path.expanduser(DEFAULT_PSP_CONFIG_PATH))
        input_gct_path = configParser.get("tests", "input_p100_subsets_gct_path")

        # Write output to the functional tests directory
        functional_tests_dir = configParser.get("tests", "functional_tests_dir")

        # Specify output file name
        OUT_NAME = "test_dry_p100_subsets_output.gct"

        args_string = ("{} -out_path {} " +
                       "-out_name {} " +
                       "-PSP_config_path {} " +
                       "-sample_nan_thresh {} " +
                       "-probe_nan_thresh {} " +
                       "-probe_sd_cutoff {} " +
                       "-dist_sd_cutoff {} " +
                       "-subset_normalize_bool " +
                       "-optim -v").format(input_gct_path, functional_tests_dir,
                                           OUT_NAME, DEFAULT_PSP_CONFIG_PATH,
                                           0.8, 0.9, 3, 3)
        args = dry.build_parser().parse_args(args_string.split())
        out_gct = dry.main(args)

        # Read in output gct created by JJ's code
        e_output_gct_path = configParser.get("tests", "output_p100_subsets_gct_path")
        e_gct = parse_gctoo.parse(e_output_gct_path)

        self.assertTrue(np.allclose(out_gct.data_df, e_gct.data_df, atol=1e-1, equal_nan=True),
                        "The data_df that was returned is wrong.")
        self.assertTrue(np.array_equal(e_gct.row_metadata_df, out_gct.row_metadata_df),
                         "The row_metadata_df that was returned is wrong.")

        # Column metadata should have one more header
        self.assertEqual(e_gct.col_metadata_df.shape[1] + 1, out_gct.col_metadata_df.shape[1],
                         ("Actual col_metadata_df should have one more header" +
                          "than e_col_metadata_df.\n" +
                          "e_gct.col_metadata_df.shape: {}, " +
                          "out_gct.col_metadata_df.shape: {}").format(e_gct.col_metadata_df.shape,
                                                                      out_gct.col_metadata_df.shape[1]))

        # Check that column metadata is correct for all but the new header
        self.assertTrue(np.array_equal(e_gct.col_metadata_df, out_gct.col_metadata_df.iloc[:, :-1]),
                        "The col_metadata_df that was returned is wrong.")

        logger.info("Mean difference between gcts: {}".format(np.nanmean(e_gct.data_df.values - out_gct.data_df.values)))
        logger.info("Mean absolute difference between gcts: {}".format(np.nanmean(np.absolute(e_gct.data_df.values - out_gct.data_df.values))))

        # Clean up
        os.remove(os.path.join(functional_tests_dir, OUT_NAME))

    def test_read_gct_and_check_provenance_code(self):
        # Get assets from config file
        configParser = ConfigParser.RawConfigParser()
        configParser.read(os.path.expanduser(DEFAULT_PSP_CONFIG_PATH))
        input_p100_gct_path = configParser.get("tests", "input_p100_gct_path")

        e_prov_code = ["PRM", "L2X"]
        e_data_df_shape = (96, 96)
        (out_gct, out_prov_code) = dry.read_gct_and_check_provenance_code(DEFAULT_PSP_CONFIG_PATH, input_p100_gct_path)
        self.assertEqual(out_prov_code, e_prov_code,
                         ("The expected provenance code is {}, " +
                          "not {}").format(e_prov_code, out_prov_code))

        self.assertEqual(out_gct.data_df.shape, e_data_df_shape,
                         ("The expected shape of the data matrix is {}, " +
                          "not {}").format(e_data_df_shape, out_gct.data_df.shape))


    def test_extract_prov_code(self):
        col_meta_df = pd.DataFrame.from_dict({"foo": ["a", "b", "c"],
                                              "provenance_code": ["PRM+L2X",
                                                                  "PRM+L2X",
                                                                  "PRM+L2X"]})
        e_prov_code = ["PRM", "L2X"]
        prov_code = dry.extract_prov_code(col_meta_df)
        self.assertEqual(e_prov_code, prov_code,
                        ("The expected provenance code is {}, " +
                         "not {}").format(e_prov_code, prov_code))

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
        e_offsets = np.array([-4.42, 2.83, 0.10], dtype=float)
        e_dists = np.array([36.62, 12.04, 0.06], dtype=float)
        (out_df, offsets, dists, success_bools) = dry.calculate_distances_and_optimize(df, (-7,7))
        self.assertTrue(np.allclose(offsets, e_offsets, atol=1e-2),
                        ("\nExpected offsets:\n{} " +
                         "\nActual offsets:\n{}").format(e_offsets, offsets))
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
        (_, _, _, _, out_prov_code) = dry.calculate_distances_and_optimize_if_needed(
            df, optim_bool, (-7, 7), prov_code)
        self.assertTrue(np.array_equal(e_prov_code, out_prov_code), (
            "The expected output provenance code is {}, not {}.".format(e_prov_code, out_prov_code)))

        # optim_bool = False, assay_type = P100
        optim_bool = False
        prov_code = ["PR1", "L2X"]
        e_prov_code = prov_code
        (_, _, _, _, out_prov_code) = dry.calculate_distances_and_optimize_if_needed(
            df, optim_bool, (-7, 7), prov_code)
        self.assertTrue(np.array_equal(e_prov_code, out_prov_code), (
            "The expected output provenance code is {}, not {}.".format(e_prov_code, out_prov_code)))

        # optim_bool = True, assay_type = GCP
        optim_bool = True
        prov_code = ["GR1", "L2X"]
        e_prov_code = prov_code
        (_, _, _, _, out_prov_code) = dry.calculate_distances_and_optimize_if_needed(
            df, optim_bool, (-7, 7), prov_code)
        self.assertTrue(np.array_equal(e_prov_code, out_prov_code), (
            "The expected output provenance code is {}, not {}.".format(e_prov_code, out_prov_code)))

        # optim_bool = True, assay_type = P100
        optim_bool = True
        prov_code = ["DIA1", "L2X"]
        e_prov_code = np.append(prov_code, "LLB")
        (_, _, _, _, out_prov_code) = dry.calculate_distances_and_optimize_if_needed(
            df, optim_bool, (-7, 7), prov_code)
        self.assertTrue(np.array_equal(e_prov_code, out_prov_code), (
            "The expected output provenance code is {}, not {}.".format(e_prov_code, out_prov_code)))



    def test_remove_sample_outliers(self):
        df = pd.DataFrame(np.array([[10, -3, 1.2, 0.6],
                                    [0.45, 0.2, 0, 0.2],
                                    [4.5, np.nan, 0.3, 0.4]], dtype=float))
        offsets = np.array([1, 2, 3, 4])
        dists = np.array([0.2, 5, 0.5, 0.4], dtype=float)
        bools = np.array([True, True, True, False], dtype=bool)
        (out, out_offsets) = dry.remove_sample_outliers(df, offsets, dists, bools,
                                         dist_sd_cutoff=1)
        e_out = df.iloc[:, [0, 2]]
        e_out_offsets = np.array([1, 3])
        self.assertTrue(out.shape == e_out.shape, (
            "expected_out.shape: {} not the same " +
            "as actual_out.shape: {}").format(e_out.shape, out.shape))
        self.assertTrue(np.array_equal(out_offsets, e_out_offsets),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_out_offsets, out_offsets))


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


    def test_make_norm_ndarray(self):
        row_df = pd.DataFrame(np.array([["8350", "1"], ["8350", "1"],
                                        ["8350", "2"], ["8350", "2"]]),
                              index=["r1", "r2", "r3", "r4"],
                              columns=["pr_gene_id", "pr_probe_normalization_group"])
        col_df = pd.DataFrame(np.array([["G-0022", "1,1"], ["G-0022", "1,1"], ["G-0022", "1,2"],
                                        ["G-0022", "2,2"], ["G-0022", "2,2"]]),
                              index=["c1", "c2", "c3", "c4", "c5"],
                              columns=["det_plate", "det_normalization_group_vector"])
        e_norm_ndarray = np.array([[1, 1, 1, 2, 2],
                                   [1, 1, 1, 2, 2],
                                   [1, 1, 2, 2, 2],
                                   [1, 1, 2, 2, 2]])
        norm_ndarray = dry.make_norm_ndarray(row_df, col_df)
        self.assertTrue(np.array_equal(norm_ndarray, e_norm_ndarray),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_norm_ndarray, norm_ndarray))

    def test_iterate_over_norm_ndarray_and_normalize(self):
        data_df = pd.DataFrame(np.array([[7, 8, 3, 8, 9],
                                         [9, 7, 4, 9, 2],
                                         [8, 6, 7, 8, 2],
                                         [4, 8, 5, 5, 7]]),
                                        index = ["a","b","c","d"],
                                        dtype="float")
        norm_ndarray = np.array([[1, 1, 1, 2, 2],
                                 [1, 1, 1, 2, 2],
                                 [1, 1, 2, 2, 2],
                                 [1, 1, 2, 2, 2]])
        e_df = pd.DataFrame(np.array([[0, 1, -4, -0.5, 0.5],
                                      [2, 0, -3, 3.5, -3.5],
                                      [1, -1, 0, 1, -5],
                                      [-2, 2, 0, 0, 2]], dtype="float"))
        out_df = dry.iterate_over_norm_ndarray_and_normalize(data_df, norm_ndarray)
        self.assertTrue(np.array_equal(out_df, e_df),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_df, out_df))

        # Slightly different norm_ndarray
        norm_ndarray = np.array([[1, 1, 1, 2, 2],
                                 [1, 1, 1, 2, 2],
                                 [1, 1, 2, 2, 3],
                                 [1, 1, 2, 2, 3]])
        e_df = pd.DataFrame(np.array([[0, 1, -4, -0.5, 0.5],
                                      [2, 0, -3, 3.5, -3.5],
                                      [1, -1, -0.5, 0.5, 0],
                                      [-2, 2, 0, 0, 0]], dtype="float"))
        out_df = dry.iterate_over_norm_ndarray_and_normalize(data_df, norm_ndarray)
        self.assertTrue(np.array_equal(out_df, e_df),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_df, out_df))

        # Totally weird but acceptable norm_ndarray
        norm_ndarray = np.array([[2, 2, 3, 3, 3],
                                 [1, 1, 2, 2, 2],
                                 [-1, -1, -1, -1, -1],
                                 [1, 1, 0, 0, 0]])
        e_df = pd.DataFrame(np.array([[-0.5, 0.5, -5, 0, 1],
                                      [1, -1, 0, 5, -2],
                                      [1, -1, 0, 1, -5],
                                      [-2, 2, 0, 0, 2]], dtype="float"))
        out_df = dry.iterate_over_norm_ndarray_and_normalize(data_df, norm_ndarray)
        self.assertTrue(np.array_equal(out_df, e_df),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_df, out_df))

    def test_subset_normalize(self):
        row_df = pd.DataFrame(np.array([["8350", "1"], ["8350", "1"],
                                        ["8350", "2"], ["8350", "2"]]),
                              index=["r1", "r2", "r3", "r4"],
                              columns=["pr_gene_id", "pr_probe_normalization_group"])
        col_df = pd.DataFrame(np.array([["G-0022", "1,1"], ["G-0022", "1,1"], ["G-0022", "1,2"],
                                        ["G-0022", "2,2"], ["G-0022", "2,2"]]),
                              index=["c1", "c2", "c3", "c4", "c5"],
                              columns=["det_plate", "det_normalization_group_vector"])
        data_df = pd.DataFrame(np.array([[7, 8, 3, 8, 9],
                                         [9, 7, 4, 9, 2],
                                         [8, 6, 7, 8, 2],
                                         [4, 8, 5, 5, 7]], dtype="float"))
        prov_code = ["GR1", "L2X", "SF3"]
        e_prov_code = ["GR1", "L2X", "SF3", "GMN"]
        e_df = pd.DataFrame(np.array([[0, 1, -4, -0.5, 0.5],
                                      [2, 0, -3, 3.5, -3.5],
                                      [1, -1, 0, 1, -5],
                                      [-2, 2, 0, 0, 2]], dtype="float"))
        (out_df, out_prov_code) = dry.subset_normalize(data_df, row_df, col_df, prov_code)
        self.assertTrue(np.array_equal(out_df, e_df),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_df, out_df))

        # Check prov code too
        self.assertEqual(out_prov_code, e_prov_code,
                         ("The expected provenance code is {}, " +
                          "not {}").format(e_prov_code, out_prov_code))

    def test_create_output_gct(self):
        row_df = pd.DataFrame(np.array([["8350", "1"], ["8350", "1"],
                                        ["8350", "2"], ["8350", "2"]]),
                              index=["r1", "r2", "r3", "r4"],
                              columns=["pr_gene_id", "pr_probe_normalization_group"])
        col_df = pd.DataFrame(np.array([["G-0021", "PRM+L2X"], ["G-0022", "PRM+L2X"], ["G-0023", "PRM+L2X"],
                                        ["G-0024", "PRM+L2X"], ["G-0025", "PRM+L2X"]]),
                              index=["c1", "c2", "c3", "c4", "c5"],
                              columns=["det_plate", "provenance_code"])
        data_df = pd.DataFrame(np.array([[1, 2, 3],
                                         [4, 5, 6]]),
                                        index=["r1", "r3"],
                                        columns=["c2", "c4", "c5"],
                                        dtype="float")
        offsets = np.array([1.1, 2.2, 3.3], dtype=float)
        prov_code = ["PRM", "L2X", "GMN", "SF3"]
        prov_code_delimiter = "+"
        out_gct = dry.create_output_gct(data_df, row_df, col_df, offsets, prov_code, prov_code_delimiter)

        e_row_df = pd.DataFrame(np.array([["8350", "1"], ["8350", "2"]]),
                              index=["r1", "r3"],
                              columns=["pr_gene_id", "pr_probe_normalization_group"])
        e_col_df = pd.DataFrame(np.array([["G-0022", "PRM+L2X+GMN+SF3", "1.1"],
                                          ["G-0024", "PRM+L2X+GMN+SF3", "2.2"],
                                          ["G-0025", "PRM+L2X+GMN+SF3", "3.3"]]),
                              index=["c2", "c4", "c5"],
                              columns=["det_plate", "provenance_code", "optimization_offset"])
        self.assertTrue(np.array_equal(out_gct.row_metadata_df, e_row_df),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_row_df, out_gct.row_metadata_df))
        self.assertTrue(np.array_equal(out_gct.col_metadata_df, e_col_df),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_col_df, out_gct.col_metadata_df))
        self.assertTrue(np.array_equal(out_gct.data_df, data_df),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(data_df, out_gct.data_df))

    def test_slice_metadata_using_already_sliced_data_df(self):
        data = pd.DataFrame([[2,3],[5,6],[11,12]],
                            index=["a","b","d"],
                            columns=["f","g"])
        row_meta = pd.DataFrame([["rm1","rm2"],["rm3","rm4"],
                                 ["rm5","rm6"],["rm7","rm8"]],
                                index=["a","b","c","d"],
                                columns=["row_field1", "row_field2"])
        col_meta = pd.DataFrame([["cm1","cm2"],["cm3","cm4"],["cm5","cm6"]],
                                index=["e","f","g"],
                                columns=["col_field1","col_field2"])
        e_row_meta = pd.DataFrame([["rm1","rm2"],["rm3","rm4"],["rm7","rm8"]],
                                index=["a","b","d"],
                                columns=["row_field1", "row_field2"])
        e_col_meta = pd.DataFrame([["cm3","cm4"],["cm5","cm6"]],
                                index=["f","g"],
                                columns=["col_field1","col_field2"])

        (out_row, out_col) = dry.slice_metadata_using_already_sliced_data_df(data,
                                                                             row_meta,
                                                                             col_meta)
        self.assertTrue(np.array_equal(out_row, e_row_meta),
                        "Row metadata dataframe is wrong: \n{}".format(out_row))
        self.assertTrue(np.array_equal(out_col, e_col_meta),
                        "Col metadata dataframe is wrong: \n{}".format(out_col))
        
    def test_check_for_subsets(self):
        row_meta = pd.DataFrame([["rm1","rm2"],["rm3","rm4"],
                                 ["rm5","rm6"],["rm7","rm8"]],
                                index=["a","b","c","d"],
                                columns=["row_field1", "row_field2"])
        col_meta = pd.DataFrame([["cm1","cm2"],["cm3","cm4"],["cm5","cm6"]],
                                index=["e","f","g"],
                                columns=["col_field1","col_field2"])
        row_field = "row_field1"
        col_field = "col_field2"
        subsets_exist = dry.check_for_subsets(row_meta, col_meta, row_field, col_field)

        self.assertTrue(subsets_exist)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
