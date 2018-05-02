import unittest
import logging
import os
import numpy as np
import pandas as pd

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.write_gct as wg
import steep

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestSteep(unittest.TestCase):

    def test_compute_similarity_bw_two_dfs(self):
        # Happy path
        df1 = pd.DataFrame.from_dict({"a": [9, 11, 2], "b": [8, 3, 1]})
        df2 = pd.DataFrame.from_dict({"c": [13, 11, 9], "d": [1, 5, 1]})

        out_df = steep.compute_similarity_bw_two_dfs(df1, df2, "spearman")

        e_df = pd.DataFrame([[0.50, 0.866], [1.00, 0.00]], index=["a", "b"], columns=["c", "d"])
        self.assertTrue(np.allclose(out_df, e_df, atol=1e-2), (
          "out_df:\n{}\ne_df:\n{}").format(out_df, e_df))

        # If dfs don't have rows in common, exception should be raised
        df3 = pd.DataFrame([[5, 8], [3, 1], [-2, 0]], index=["a", "b", "c"])

        with self.assertRaises(AssertionError) as e:
            steep.compute_similarity_bw_two_dfs(df1, df3, "spearman")
        self.assertIn("Are you sure", str(e.exception))


    def test_compute_similarity_within_df(self):
        ### Check that NaNs are handled correctly
        METHOD = "spearman"
        df = pd.DataFrame.from_dict({"a": [3, 5, 1, 2, 4],
                                     "b": [2, 9, np.nan, 8, 1],
                                     "c": [np.nan, 3, 5, 2, 8]})
        out_df = steep.compute_similarity_within_df(df, METHOD)

        # Expect col1 v. col2 to be [3,5,2,4] v. [2,9,8,1]
        col1_v_col2_cor = pd.Series([3, 5, 2, 4]).corr(pd.Series([2, 9, 8, 1]), method=METHOD)
        self.assertAlmostEquals(col1_v_col2_cor, out_df.iloc[0, 1])

        # Expect col1 v. col3 to be [5,1,2,4] v. [3,5,2,8]
        col1_v_col3_cor = pd.Series([5, 1, 2, 4]).corr(pd.Series([3, 5, 2, 8]), method=METHOD)
        self.assertAlmostEquals(col1_v_col3_cor, out_df.iloc[0, 2])

        # Expect col2 v. col3 to be [9,8,1] v. [3,2,8]
        col2_v_col3_cor = pd.Series([9, 8, 1]).corr(pd.Series([3, 2, 8]), method=METHOD)
        self.assertAlmostEquals(col2_v_col3_cor, out_df.iloc[1, 2])

        e_df = pd.DataFrame([[1.0, 0.2, 0.0], [0.2, 1.0, -0.5], [0.0, -0.5, 1.0]],
                            index=["a", "b", "c"], columns=["a", "b", "c"])
        self.assertTrue(out_df.equals(e_df))

        ### Check error handling for bad similarity_metric
        with self.assertRaises(Exception) as e:
            steep.compute_similarity_within_df(df, "wtcs")
        self.assertIn("similarity metric must be", str(e.exception))


    def test_main(self):
        # Prepare the data
        col_meta_df_index1 = pd.Index(["a", "b"], name="cid")
        col_meta_df_index2 = pd.Index(["c", "d", "e"], name="cid")
        col_meta_df_columns = pd.Index(["chd1"], name="chd")
        col_meta_df1 = pd.DataFrame([["a1"], ["b1"]],
            index=col_meta_df_index1, columns=col_meta_df_columns)
        col_meta_df2 = pd.DataFrame([["c1"], ["d1"], ["e1"]],
            index=col_meta_df_index2, columns=col_meta_df_columns)
        row_meta_df_index = pd.Index(["A", "B", "C"], name="rid")
        row_meta_df_columns = pd.Index(["rhd1"], name="rhd")
        row_meta_df = pd.DataFrame([["rhdA"], ["rhdB"], ["rhdC"]],
            index=row_meta_df_index, columns=row_meta_df_columns)
        df1 = pd.DataFrame.from_dict({"a":[1, 3, 7], "b":[5, 2, 7]})
        df1.index = row_meta_df_index
        df1.columns.name = "rhd"
        df2 = pd.DataFrame.from_dict({"c":[9, 11, 2], "d":[4, 10, 3], "e": [7, 4, 3]})
        df2.index = row_meta_df_index
        df2.columns.name = "rhd"
        gct1 = GCToo.GCToo(df1, row_metadata_df=row_meta_df, col_metadata_df=col_meta_df1)
        gct2 = GCToo.GCToo(df2, row_metadata_df=row_meta_df, col_metadata_df=col_meta_df2)

        # What to name files
        gct1_name = "test_steep_main_gct1.gct"
        gct2_name = "test_steep_main_gct2.gct"
        out_name1 = "test_steep_main_sim1.gct"
        out_name2 = "test_steep_main_sim2.gct"
        files = [gct1_name, gct2_name, out_name1, out_name2]

        # Write the gcts
        wg.write(gct1, gct1_name)
        wg.write(gct2, gct2_name)

        # Compute similarity within gct
        arg_string1 = "-i {} -o {}".format(gct1_name, out_name1)
        args1 = steep.build_parser().parse_args(arg_string1.split())
        steep.main(args1)

        sim_gct1 = parse.parse(out_name1)
        e_index1 = pd.Index(["a", "b"])
        spearman_series = pd.Series(["spearman", "spearman"],
                                    index=["a", "b"],
                                    name="similarity_metric")
        e_col_meta_df = pd.concat([col_meta_df1, spearman_series], axis=1)
        self.assertTrue(sim_gct1.data_df.index.equals(e_index1), (
            "\nsim_gct1.data_df.index: {}\ne_index1: {}").format(
                sim_gct1.data_df.index, e_index1))
        self.assertTrue(sim_gct1.data_df.columns.equals(e_index1), (
            "\nsim_gct1.data_df.columns: {}\ne_index1: {}").format(
                sim_gct1.data_df.columns, e_index1))
        self.assertTrue(sim_gct1.row_metadata_df.equals(e_col_meta_df), (
            "\nsim_gct1.row_metadata_df:\n{}\ne_col_meta_df:\n{}").format(
                sim_gct1.row_metadata_df, e_col_meta_df))
        self.assertTrue(sim_gct1.col_metadata_df.equals(e_col_meta_df), (
            "\nsim_gct1.col_metadata_df:\n{}\ne_col_meta_df:\n{}").format(
                sim_gct1.col_metadata_df, e_col_meta_df))

        # Compute similarity across gcts
        arg_string2 = "-i {} -o {} -i2 {}".format(gct1_name,
            out_name2, gct2_name)
        args2 = steep.build_parser().parse_args(arg_string2.split())
        steep.main(args2)

        sim_gct2 = parse.parse(out_name2)
        e_index2 = pd.Index(["a", "b"])
        e_columns2 = pd.Index(["c", "d", "e"])
        spearman_series2 = pd.Series(["spearman", "spearman", "spearman"],
                                    index=["c", "d", "e"],
                                    name="similarity_metric")
        e_col_meta_df2 = pd.concat([col_meta_df2, spearman_series2], axis=1)
        self.assertTrue(sim_gct2.data_df.index.equals(e_index2), (
            "\nsim_gct2.data_df.index: {}\ne_index2: {}").format(
                sim_gct2.data_df.index, e_index2))
        self.assertTrue(sim_gct2.data_df.columns.equals(e_columns2), (
            "\nsim_gct2.data_df.columns: {}\ne_columns2: {}").format(
                sim_gct2.data_df.columns, e_columns2))
        self.assertTrue(sim_gct2.row_metadata_df.equals(e_col_meta_df), (
            "\nsim_gct2.row_metadata_df:\n{}\ne_col_meta_df:\n{}").format(
                sim_gct2.row_metadata_df, e_col_meta_df))
        self.assertTrue(sim_gct2.col_metadata_df.equals(e_col_meta_df2), (
            "\nsim_gct2.col_metadata_df:\n{}\ne_col_meta_df2:\n{}").format(
                sim_gct2.col_metadata_df, e_col_meta_df2))

        # Remove the gcts and sim outputs
        for file in files:
            os.remove(file)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
