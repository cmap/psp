import unittest
import logging
import setup_GCToo_logger as setup_logger
import ConfigParser

import os
import pandas as pd
import numpy as np
import parse_gctoo as pg

logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Read config file
configParser = ConfigParser.RawConfigParser()
config_path = os.path.expanduser("~/.PSP_config")
configParser.read(config_path)

# Get location of functional test files
functional_tests_path = configParser.get("tests", "functional_tests_path")
gct_filename = "LJP.gct"


"""
TO-DO(lev): use custom dataframes instead???
"""


class TestParseGCToo(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.gct_filepath = "/Users/lev/code/test_LJP.gct"
        self.gct_dims = [978, 377, 11, 35]
        self.full_df = pd.DataFrame([["chd1", "", "", "a", "b"],
                                          ["chd2", "", "", "55", "61"],
                                          ["chd3", "", "", "nah", "nope"],
                                          ["rid1", "C", "D", "0.3", "0.2"],
                                          ["rid2", "1.0", "2.0", np.nan, "0.9"]],
                                         columns=["id", "rhd1", "rhd2", "cid1", "cid2"])
        self.full_df_dims = [2, 2, 2, 3]

    def test_read_version_and_dims(self):
        version = "#1.3"
        dims = ["10", "15", "3", "4"]
        fname = "testing_testing"

        f = open(fname, "wb")
        f.write((version + "\n"))
        f.write((dims[0] + "\t" + dims[1] + "\t" + dims[2] + "\t" + dims[3] + "\n"))
        f.close()

        (actual_version, n_rows, n_cols, n_metarows, n_metacols) = pg.read_version_and_dims(fname)
        self.assertEqual(actual_version, version)
        self.assertEqual(n_rows, int(dims[0]))
        self.assertEqual(n_metacols, int(dims[3]))

        # Remove the file I created
        os.remove(fname)

    def test_parse_into_3_df(self):
        dims = self.gct_dims
        (row_metadata, col_metadata, data) = pg.parse_into_3_df(self.gct_filepath,
                                                                dims[0], dims[1], dims[2],
                                                                dims[3], None)

        # Check shapes of outputs
        self.assertTrue(row_metadata.shape == (dims[0], dims[2]),
                        ("row_metadata.shape = {} " +
                         "but expected it to be ({}, {})").format(row_metadata.shape,
                                                                  dims[0], dims[2]))
        self.assertTrue(col_metadata.shape == (dims[3], dims[1]),
                        ("row_metadata.shape = {} " +
                         "but expected it to be ({}, {})").format(row_metadata.shape,
                                                                  dims[3], dims[1]))
        self.assertTrue(data.shape == (dims[0], dims[1]),
                        ("row_metadata.shape = {} " +
                         "but expected it to be ({}, {})").format(row_metadata.shape,
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
        self.assertTrue(col_metadata.iloc[0, dims[1] - 1] == correct_str,
                        ("The last value in the first row of column metadata should be " +
                         correct_str + " not {}").format(col_metadata.iloc[0, dims[1] - 1]))

        # Check headers
        correct_str = "LJP005_A375_24H_X1_B19:P24"
        self.assertTrue(list(col_metadata)[dims[1] - 1] == correct_str,
                        ("The last column metadata header should be " +
                         correct_str + " not {}").format(list(col_metadata)[dims[1] - 1]))
        correct_str = "bead_batch"
        self.assertTrue(col_metadata.index.values[3] == correct_str,
                        ("The fourth column metadata index value should be " +
                         correct_str + " not {}").format(col_metadata.index.values[3]))
        correct_str = "203897_at"
        self.assertTrue(row_metadata.index.values[dims[0] - 1] == correct_str,
                        ("The last row metadata index value should be " + correct_str +
                         " not {}").format(row_metadata.index.values[dims[0] - 1]))
        self.assertTrue(data.index.values[dims[0] - 1] == correct_str,
                        ("The last data index value should be " + correct_str +
                         " not {}").format(data.index.values[dims[0] - 1]))

    # def test_assemble_row_metadata(self):
    #     pg.assemble_row_metadata(full_df, num_col_metadata, num_data_rows, num_row_metadata)
    #
    # def test_assemble_col_metadata(self):
    #     pg.assemble_col_metadata(full_df, num_col_metadata, num_row_metadata, num_data_cols)
    #
    def test_assemble_data(self):
        dims = self.full_df_dims
        data = pg.assemble_data(self.full_df, dims[3], dims[0], dims[2], dims[1])
        self.assertTrue(np.isnan(data.iloc[1, 0]),
                        ("The 2nd row, 1st column data point should be nan, " +
                         "not {}").format(data.iloc[1, 0]))
    #
    # def test_create_gctoo_obj(self):
    #     pg.test_create_gctoo_obj(gct_filename, version, row_metadata_df, col_metadata_df, data_df)
    #
    # def test_main(self):
    #     pg.main(gct_filename)


if __name__ == "__main__":
    # setup_logger.setup(verbose=True)
    unittest.main()
