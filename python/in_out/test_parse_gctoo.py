import unittest
import logging
import setup_GCToo_logger as setup_logger
import ConfigParser

import os
import pandas as pd
import numpy as np
import parse_gctoo as pg

"""
N.B. To run the functional tests, you will have to create a config file
called .gctoo_config in your home directory.
"""

# Set up logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Config file location
DEFAULT_CONFIG_PATH = "~/.gctoo_config"

# Read config file if it exists
config_path = os.path.expanduser(DEFAULT_CONFIG_PATH)
if os.path.exists(config_path):
    configParser = ConfigParser.RawConfigParser()
    configParser.read(config_path)
    functional_tests_path = configParser.get("tests", "functional_tests_path")

# Otherwise, skip functional tests
else:
    logger.info(("No config file could be found at config_path: {}. " +
                 "No functional tests will be run.").format(config_path))
    functional_tests_path = None


class TestParseGCToo(unittest.TestCase):
    def test_read_version_and_dims(self):
        version = "1.3"
        dims = ["10", "15", "3", "4"]
        fname = "testing_testing"

        f = open(fname, "wb")
        f.write(("#" + version + "\n"))
        f.write((dims[0] + "\t" + dims[1] + "\t" + dims[2] + "\t" + dims[3] + "\n"))
        f.close()

        (actual_version, n_rows, n_cols, n_rhd, n_chd) = pg.read_version_and_dims(fname)
        self.assertEqual(actual_version, version)
        self.assertEqual(n_rows, int(dims[0]))
        self.assertEqual(n_chd, int(dims[3]))

        # Remove the file I created
        os.remove(fname)

    # Only run test if config file exists
    @unittest.skipUnless(functional_tests_path is not None,
                         "Config file must exist to run functional tests.")
    def test_parse_into_3_df(self):
        # This test uses real data
        gct_filepath = os.path.join(functional_tests_path, "LJP.gct")
        dims = [978, 377, 11, 35]
        (row_metadata, col_metadata, data) = pg.parse_into_3_df(gct_filepath,
                                                                dims[0], dims[1], dims[2],
                                                                dims[3], None)

        # Check shapes of outputs
        self.assertTrue(row_metadata.shape == (dims[0], dims[2]),
                        ("row_metadata.shape = {} " +
                         "but expected it to be ({}, {})").format(row_metadata.shape,
                                                                  dims[0], dims[2]))
        self.assertTrue(col_metadata.shape == (dims[1], dims[3]),
                        ("col_metadata.shape = {} " +
                         "but expected it to be ({}, {})").format(col_metadata.shape,
                                                                  dims[1], dims[3]))
        self.assertTrue(data.shape == (dims[0], dims[1]),
                        ("data.shape = {} " +
                         "but expected it to be ({}, {})").format(data.shape,
                                                                  dims[0], dims[1]))

        # Type-check the data
        self.assertTrue(isinstance(data.iloc[0, 0], float), "The data should be a float, not a string.")

        # Check a few values
        correct_val = 11.3819
        self.assertTrue(data.iloc[0, 0] == correct_val,
                        ("The first value in the data matrix should be " +
                         str(correct_val) + " not {}").format(data.iloc[0, 0]))
        correct_val = 5.1256
        self.assertTrue(data.iloc[dims[0] - 1, dims[1] - 1] == correct_val,
                        ("The last value in the data matrix should be " +
                         str(correct_val) + " not {}").format(data.iloc[dims[0] - 1, dims[1] - 1]))
        correct_str = "LUA-4000"
        self.assertTrue(row_metadata.iloc[2, 3] == correct_str,
                        ("The 3rd row, 4th column of the row metadata should be " +
                         correct_str + " not {}").format(row_metadata.iloc[2, 3]))
        correct_str = "57"
        self.assertTrue(col_metadata.iloc[dims[1] - 1, 0] == correct_str,
                        ("The last value in the first column of column metadata should be " +
                         correct_str + " not {}").format(col_metadata.iloc[dims[1] - 1, 0]))

        # Check headers
        correct_str = "LJP005_A375_24H_X1_B19:P24"
        self.assertTrue(col_metadata.index.values[dims[1] - 1] == correct_str,
                        ("The last column metadata header should be " +
                         correct_str + " not {}").format(col_metadata.index.values[dims[1] - 1]))
        correct_str = "bead_batch"
        self.assertTrue(list(col_metadata)[3] == correct_str,
                        ("The fourth column metadata index value should be " +
                         correct_str + " not {}").format(list(col_metadata)[3]))
        correct_str = "203897_at"
        self.assertTrue(row_metadata.index.values[dims[0] - 1] == correct_str,
                        ("The last row metadata index value should be " + correct_str +
                         " not {}").format(row_metadata.index.values[dims[0] - 1]))
        self.assertTrue(data.index.values[dims[0] - 1] == correct_str,
                        ("The last data index value should be " + correct_str +
                         " not {}").format(data.index.values[dims[0] - 1]))

    def test_assemble_row_metadata(self):
        full_df = pd.DataFrame(
            [["chd1", "", "", "a", "b"],
             ["chd2", "", "", "55", "61"],
             ["chd3", "", "", "nah", "nope"],
             ["rid1", "C", "D", "0.3", "0.2"],
             ["rid2", "1.0", "2.0", np.nan, "0.9"]],
            columns=["id", "rhd1", "rhd2", "cid1", "cid2"])
        full_df_dims = [2, 2, 2, 3]
        e_row_df = pd.DataFrame([["C", "D"], ["1.0", "2.0"]],
                                index=["rid1", "rid2"],
                                columns=["rhd1", "rhd2"])
        row_df = pg.assemble_row_metadata(full_df, full_df_dims[3],
                                          full_df_dims[0], full_df_dims[2])
        self.assertTrue(row_df.equals(e_row_df))

    def test_assemble_col_metadata(self):
        full_df = pd.DataFrame(
            [["chd1", "", "", "a", "b"],
             ["chd2", "", "", "55", "61"],
             ["chd3", "", "", "nah", "nope"],
             ["rid1", "C", "D", "0.3", "0.2"],
             ["rid2", "1.0", "2.0", np.nan, "0.9"]],
            columns=["id", "rhd1", "rhd2", "cid1", "cid2"])
        full_df_dims = [2, 2, 2, 3]
        e_col_df = pd.DataFrame([["a", "55", "nah"], ["b", "61", "nope"]],
                                index=["cid1", "cid2"],
                                columns=["chd1", "chd2", "chd3"])
        col_df = pg.assemble_col_metadata(full_df, full_df_dims[3],
                                          full_df_dims[2], full_df_dims[1])
        self.assertTrue(col_df.equals(e_col_df))

    def test_assemble_data(self):
        full_df = pd.DataFrame(
            [["chd1", "", "", "a", "b"],
             ["chd2", "", "", "55", "61"],
             ["chd3", "", "", "nah", "nope"],
             ["rid1", "C", "D", "0.3", "0.2"],
             ["rid2", "1.0", "2.0", np.nan, "0.9"]],
            columns=["id", "rhd1", "rhd2", "cid1", "cid2"])
        full_df_dims = [2, 2, 2, 3]
        e_data_df = pd.DataFrame([[0.3, 0.2], [np.nan, 0.9]],
                                 index=["rid1", "rid2"],
                                 columns=["cid1", "cid2"])
        data_df = pg.assemble_data(full_df, full_df_dims[3], full_df_dims[0],
                                   full_df_dims[2], full_df_dims[1])
        self.assertTrue(data_df.equals(e_data_df))

    # Only run test if config file exists
    @unittest.skipUnless(functional_tests_path is not None,
                         "Config file must exist to run functional tests.")
    def test_parse(self):
        # L1000 gct
        LJP_file_path = os.path.join(functional_tests_path, "LJP.gct")
        LJP_gct = pg.parse(LJP_file_path)

        # Check a few values
        correct_val = 11.3819
        self.assertTrue(LJP_gct.data_df.iloc[0, 0] == correct_val,
                        ("The first value in the data matrix should be " +
                         "{} not {}").format(str(correct_val), LJP_gct.data_df.iloc[0, 0]))
        correct_val = "58"
        self.assertTrue(LJP_gct.col_metadata_df.iloc[0, 0] == correct_val,
                        ("The first value in the column metadata should be " +
                         "{} not {}").format(str(correct_val), LJP_gct.col_metadata_df.iloc[0, 0]))
        correct_val = "Analyte 11"
        self.assertTrue(LJP_gct.row_metadata_df.iloc[0, 0] == correct_val,
                        ("The first value in the row metadata should be " +
                         "{} not {}").format(str(correct_val), LJP_gct.row_metadata_df.iloc[0, 0]))

        # PSP gct
        p100_file_path = os.path.join(functional_tests_path, "p100_prm_plate29_3H.gct")
        p100_gct = pg.parse(p100_file_path)

        # Check a few values
        correct_val = 0.918157217057044
        self.assertTrue(p100_gct.data_df.iloc[0, 0] == correct_val,
                        ("The first value in the data matrix should be " +
                         "{} not {}").format(str(correct_val), p100_gct.data_df.iloc[0, 0]))
        correct_val = "MCF7"
        self.assertTrue(p100_gct.col_metadata_df.iloc[0, 0] == correct_val,
                        ("The first value in the column metadata should be " +
                         "{} not {}").format(str(correct_val), p100_gct.col_metadata_df.iloc[0, 0]))
        correct_val = "1859"
        self.assertTrue(p100_gct.row_metadata_df.iloc[0, 0] == correct_val,
                        ("The first value in the row metadata should be " +
                         "{} not {}").format(str(correct_val), p100_gct.row_metadata_df.iloc[0, 0]))


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
