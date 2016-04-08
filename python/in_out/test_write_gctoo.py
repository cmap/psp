import logging
import setup_GCToo_logger as setup_logger
import unittest
import ConfigParser
import os
import pandas as pd
import parse_gctoo
import write_gctoo

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


class TestWriteGCToo(unittest.TestCase):
    def test_write_version_and_dims(self):
        # Write
        fname = "test_file.gct"
        f = open(fname, "wb")
        write_gctoo.write_version_and_dims("1.3", ["1", "2", "3", "4"], f)
        f.close()

        # Read and then remove
        f = open(fname, "r")
        version_string = f.readline().strip()
        dims = f.readline().strip().split("\t")
        f.close()
        os.remove(fname)

        # Check that it was written correctly
        self.assertEqual(version_string, "#1.3")
        self.assertEqual(dims, ["1", "2", "3", "4"])

    def test_assemble_full_df(self):
        row_meta_df = pd.DataFrame([["Analyte 11", "11", "dp52"],
                                    ["Analyte 12", "12", "dp52"]],
                                   index=["200814_at", "218597_s_at"],
                                   columns=["pr_analyte_id", "pr_analyte_num", "pr_bset_id"])
        col_meta_df = pd.DataFrame([["8.38", "-666", "DMSO", "24 h"],
                                    ["7.7", "-666", "DMSO", "24 h"],
                                    ["8.18", "-666", "DMSO", "24 h"]],
                                   index=["LJP005_A375_24H_X1_B19:A03",
                                          "LJP005_A375_24H_X1_B19:A04",
                                          "LJP005_A375_24H_X1_B19:A05"],
                                   columns=["qc_iqr", "pert_idose", "pert_iname", "pert_itime"])
        data_df = pd.DataFrame([[11.3819, 11.3336, 11.4486],
                                [10.445, 10.445, 10.3658]],
                               index=["200814_at", "218597_s_at"],
                               columns=["LJP005_A375_24H_X1_B19:A03",
                                        "LJP005_A375_24H_X1_B19:A04",
                                        "LJP005_A375_24H_X1_B19:A05"])
        e_df = pd.DataFrame(
            [['qc_iqr', '-666', '-666', '-666', '8.38', '7.7', '8.18'],
             ['pert_idose', '-666', '-666', '-666', '-666', '-666', '-666'],
             ['pert_iname', '-666', '-666', '-666', 'DMSO', 'DMSO', 'DMSO'],
             ['pert_itime', '-666', '-666', '-666', '24 h', '24 h', '24 h'],
             ['200814_at', 'Analyte 11', '11', 'dp52', 11.3819, 11.3336, 11.4486],
             ['218597_s_at', 'Analyte 12', '12', 'dp52', 10.445, 10.445, 10.3658]],
            columns=['id', 'pr_analyte_id', 'pr_analyte_num', 'pr_bset_id',
                     'LJP005_A375_24H_X1_B19:A03', 'LJP005_A375_24H_X1_B19:A04',
                     'LJP005_A375_24H_X1_B19:A05'])
        full_df = write_gctoo.assemble_full_df(row_meta_df, col_meta_df, data_df, "NaN", "-666")

        self.assertTrue(full_df.equals(e_df))
        self.assertEqual(full_df.columns.values[0], "id")
        self.assertEqual(full_df.iloc[4, 0], "200814_at")
        self.assertEqual(full_df.ix[2, "pr_bset_id"], "-666")

    def test_write_full_df(self):
        full_df = pd.DataFrame(
            [['qc_iqr', '-666', '-666', '-666', '8.38', '7.7', '8.18'],
             ['pert_idose', '-666', '-666', '-666', '-666', '-666', '-666'],
             ['pert_iname', '-666', '-666', '-666', 'DMSO', 'DMSO', 'DMSO'],
             ['pert_itime', '-666', '-666', '-666', '24 h', '24 h', '24 h'],
             ['200814_at', 'Analyte 11', '11', 'dp52', 11.3819, 11.3336, 11.4486],
             ['218597_s_at', 'Analyte 12', '12', 'dp52', 9.5063, 10.445, 10.3658]],
            columns=['id', 'pr_analyte_id', 'pr_analyte_num', 'pr_bset_id',
                     'LJP005_A375_24H_X1_B19:A03', 'LJP005_A375_24H_X1_B19:A04',
                     'LJP005_A375_24H_X1_B19:A05'])

        fname = "test_file.gct"
        f = open(fname, "wb")
        write_gctoo.write_full_df(full_df, f, "NaN", None)
        f.close()
        os.remove(fname)

        f2 = open(fname, "wb")
        f2.write("#1.3\n")
        f2.write("1\t2\t3\t4\n")
        write_gctoo.write_full_df(full_df, f2, "NaN", None)
        f2.close()
        os.remove(fname)

    def test_append_dims_and_file_extension(self):
        data_df = pd.DataFrame([[1, 2], [3, 4]])
        fname_no_gct = "a/b/file"
        fname_gct = "a/b/file.gct"
        e_fname = "a/b/file_n2x2.gct"

        fname_out = write_gctoo.append_dims_and_file_extension(fname_no_gct, data_df)
        self.assertEqual(fname_out, e_fname)

        fname_out = write_gctoo.append_dims_and_file_extension(fname_gct, data_df)
        self.assertEqual(fname_out, e_fname)

    # Only run test if config file exists
    @unittest.skipUnless(functional_tests_path is not None,
                         "Config file must exist to run functional tests.")
    def test_LJP_functional(self):
        LJP_in_path = os.path.join(functional_tests_path, "LJP.gct")
        LJP_out_path = os.path.join(functional_tests_path, "test_LJP_writing.gct")

        # Read in original gct file
        LJP_in_gct = parse_gctoo.parse(LJP_in_path)

        # Read in new gct file
        write_gctoo.write(LJP_in_gct, LJP_out_path)
        LJP_out_gct = parse_gctoo.parse(LJP_out_path)

        self.assertTrue(LJP_in_gct.data_df.equals(LJP_out_gct.data_df))
        self.assertTrue(LJP_in_gct.row_metadata_df.equals(LJP_out_gct.row_metadata_df))
        self.assertTrue(LJP_in_gct.col_metadata_df.equals(LJP_out_gct.col_metadata_df))

        # Clean up
        os.remove(LJP_out_path)
        logger.info("{} has been removed.".format(LJP_out_path))

    # Only run test if config file exists
    @unittest.skipUnless(functional_tests_path is not None,
                         "Config file must exist to run functional tests.")
    def test_p100_functional(self):
        p100_in_path = os.path.join(functional_tests_path, "p100_prm_plate29_3H.gct")
        p100_out_path = os.path.join(functional_tests_path, "test_p100_writing.gct")

        # Read in original gct file
        p100_in_gct = parse_gctoo.parse(p100_in_path)

        # Read in new gct file
        write_gctoo.write(p100_in_gct, p100_out_path)
        p100_out_gct = parse_gctoo.parse(p100_out_path)

        self.assertTrue(p100_in_gct.data_df.equals(p100_out_gct.data_df))
        self.assertTrue(p100_in_gct.row_metadata_df.equals(p100_out_gct.row_metadata_df))
        self.assertTrue(p100_in_gct.col_metadata_df.equals(p100_out_gct.col_metadata_df))

        # Clean up
        os.remove(p100_out_path)
        logger.info("{} has been removed.".format(p100_out_path))

if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
