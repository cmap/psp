import unittest
import logging
import utils.setup_logger as setup_logger
import os
import numpy as np
import pandas as pd

import in_out.parse_gctoo as parse_gctoo
import in_out.GCToo as GCToo
import dry

"""
This code should be run from broadinstitute.psp.

The dry directory contains a directory called functional_tests
that has the assets required for the 3 functional tests below. For
functional tests, I just check that they run to completion.

"""
# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Functional tests dir lives within the dry directory
FUNCTIONAL_TESTS_DIR = "dry/functional_tests"

# Set to false if you want to see what output is created
CLEANUP = True

# N.B. e_out is expected_output.
class TestDryAgainstRCode(unittest.TestCase):

    def test_p100_main(self):
        INPUT_GCT_PATH = os.path.join(FUNCTIONAL_TESTS_DIR, "p100_prm_plate29_3H.gct")
        JJ_OUTPUT_GCT = os.path.join(FUNCTIONAL_TESTS_DIR, "JJ_p100_prm_plate29_3H_processed.gct")
        OUT_NAME = "test_dry_p100_output.gct"
        OUT_PW_NAME = "test_dry_p100_output.pw"

        args_string = ("{} {} " +
                       "-out_name {} " +
                       "-out_pw_name {} " +
                       "-sample_nan_thresh {} " +
                       "-probe_nan_thresh {} " +
                       "-dist_sd_cutoff {} " +
                       "-v").format(INPUT_GCT_PATH, FUNCTIONAL_TESTS_DIR,
                                    OUT_NAME, OUT_PW_NAME, 0.8, 0.9, 3)
        args = dry.build_parser().parse_args(args_string.split())
        out_gct = dry.main(args)

        # Compare to output gct created by JJ's code
        e_gct = parse_gctoo.parse(JJ_OUTPUT_GCT)

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

        # Check that column metadata is correct (ignore last 3 entries so that we can skip provenance code)
        self.assertTrue(np.array_equal(e_gct.col_metadata_df.iloc[:, :-2], out_gct.col_metadata_df.iloc[:, :-3]),
                        "The col_metadata_df that was returned is wrong.")

        logger.info("Mean difference between gcts: {}".format(np.nanmean(e_gct.data_df.values - out_gct.data_df.values)))
        logger.info("Mean absolute difference between gcts: {}".format(np.nanmean(np.absolute(e_gct.data_df.values - out_gct.data_df.values))))

        # Clean up
        if CLEANUP:
            os.remove(os.path.join(FUNCTIONAL_TESTS_DIR, OUT_NAME))
            os.remove(os.path.join(FUNCTIONAL_TESTS_DIR, OUT_PW_NAME))

    def test_GCP_main(self):
        INPUT_GCT_PATH = os.path.join(FUNCTIONAL_TESTS_DIR, "gcp_gr1_plate31.gct")
        JJ_OUTPUT_GCT = os.path.join(FUNCTIONAL_TESTS_DIR, "JJ_gcp_gr1_plate31_processed.gct")
        OUT_NAME = "test_dry_gcp_output.gct"
        OUT_PW_NAME = "test_dry_gcp_output.pw"

        args_string = ("{} {} " +
                       "-out_name {} " +
                       "-out_pw_name {} " +
                       "-sample_nan_thresh {} " +
                       "-probe_nan_thresh {} " +
                       "-probe_sd_cutoff {} " +
                       "-v").format(INPUT_GCT_PATH, FUNCTIONAL_TESTS_DIR,
                                    OUT_NAME, OUT_PW_NAME, 0.5, 0.5, 4)
        args = dry.build_parser().parse_args(args_string.split())
        out_gct = dry.main(args)

        # Compare to output gct created by JJ's code
        e_gct = parse_gctoo.parse(JJ_OUTPUT_GCT)

        self.assertTrue(np.allclose(out_gct.data_df, e_gct.data_df, atol=1e-1, equal_nan=True),
                        "The data_df that was returned is wrong.")
        self.assertTrue(np.array_equal(e_gct.row_metadata_df, out_gct.row_metadata_df),
                        "The row_metadata_df that was returned is wrong.")

        # Check that column metadata is correct (ignore last 2 entries so that we can skip provenance code)
        self.assertTrue(np.array_equal(e_gct.col_metadata_df.iloc[:, :-2], out_gct.col_metadata_df.iloc[:, :-2]),
                        "The col_metadata_df that was returned is wrong.")

        logger.info("Mean difference between gcts: {}".format(np.nanmean(e_gct.data_df.values - out_gct.data_df.values)))
        logger.info("Mean absolute difference between gcts: {}".format(np.nanmean(np.absolute(e_gct.data_df.values - out_gct.data_df.values))))

        # No samples should have been filtered out
        out_pw_df = pd.read_csv(
            os.path.join(FUNCTIONAL_TESTS_DIR, OUT_PW_NAME), sep="\t")
        self.assertTrue(all(out_pw_df["remains_after_poor_coverage_filtration"]))

        # Clean up
        if CLEANUP:
            os.remove(os.path.join(FUNCTIONAL_TESTS_DIR, OUT_NAME))
            os.remove(os.path.join(FUNCTIONAL_TESTS_DIR, OUT_PW_NAME))


    def test_p100_subset_main(self):
        INPUT_GCT_PATH = os.path.join(FUNCTIONAL_TESTS_DIR, "p100_prm_plate35_subsets.gct")
        JJ_OUTPUT_GCT = os.path.join(FUNCTIONAL_TESTS_DIR, "JJ_p100_prm_plate35_subsets_processed.gct")
        OUT_NAME = "test_dry_p100_subsets_output.gct"
        OUT_PW_NAME = "test_dry_p100_subsets_output.pw"

        args_string = ("{} {} " +
                       "-out_name {} " +
                       "-out_pw_name {} " +
                       "-sample_nan_thresh {} " +
                       "-probe_nan_thresh {} " +
                       "-probe_sd_cutoff {} " +
                       "-dist_sd_cutoff {} " +
                       "-v").format(INPUT_GCT_PATH, FUNCTIONAL_TESTS_DIR,
                                    OUT_NAME, OUT_PW_NAME, 0.8, 0.9, 3, 3)
        args = dry.build_parser().parse_args(args_string.split())
        out_gct = dry.main(args)

        # Compare to output gct created by JJ's code
        e_gct = parse_gctoo.parse(JJ_OUTPUT_GCT)

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

        # Check that column metadata is correct for all but the new header and provenance code
        self.assertTrue(np.array_equal(e_gct.col_metadata_df.iloc[:, :-1], out_gct.col_metadata_df.iloc[:, :-2]),
                        "The col_metadata_df that was returned is wrong.")

        logger.info("Mean difference between gcts: {}".format(np.nanmean(e_gct.data_df.values - out_gct.data_df.values)))
        logger.info("Mean absolute difference between gcts: {}".format(np.nanmean(np.absolute(e_gct.data_df.values - out_gct.data_df.values))))

        # Clean up
        if CLEANUP:
            os.remove(os.path.join(FUNCTIONAL_TESTS_DIR, OUT_NAME))
            os.remove(os.path.join(FUNCTIONAL_TESTS_DIR, OUT_PW_NAME))


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
