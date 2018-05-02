import unittest
import logging
import os
import numpy as np
import pandas as pd
import scipy.stats as stats

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.GCToo as GCToo
import sip

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

FUNCTIONAL_TESTS_DIR = "sip/functional_tests"


class TestSip(unittest.TestCase):

    def test_main(self):
        test_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_sip_in_test.gct")
        bg_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_sip_in_bg.gct")
        out_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_sip_main_out.gct")

        args_string = "-t {} -b {} -o {} -tfq {} -tft {} -bf {} -s {}".format(
            test_gct_path, bg_gct_path, out_path, "pert_iname",
            "pert_iname", "pert_iname", "|")
        args = sip.build_parser().parse_args(args_string.split())

        # Run main method
        sip.main(args)

        # Compare the output of main with the expected output
        e_out_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_sip_expected_conn.gct")
        e_out_gct = parse.parse(e_out_path)
        out_gct = parse.parse(out_path)

        logger.debug("e_out_gct.data_df:\n{}".format(e_out_gct.data_df))
        logger.debug("out_gct.data_df:\n{}".format(out_gct.data_df))
        pd.util.testing.assert_frame_equal(e_out_gct.data_df, out_gct.data_df,
                                           check_less_precise=3)

        logger.debug("e_out_gct.row_metadata_df:\n{}".format(e_out_gct.row_metadata_df))
        logger.debug("out_gct.row_metadata_df:\n{}".format(out_gct.row_metadata_df))
        pd.util.testing.assert_frame_equal(
            e_out_gct.row_metadata_df, out_gct.row_metadata_df)

        logger.debug("e_out_gct.col_metadata_df:\n{}".format(e_out_gct.col_metadata_df))
        logger.debug("out_gct.col_metadata_df:\n{}".format(out_gct.col_metadata_df))
        pd.util.testing.assert_frame_equal(
            e_out_gct.col_metadata_df, out_gct.col_metadata_df)

        # Remove the created file
        os.remove(out_path)

    def test_check_symmetry(self):
        df_mat = np.random.randn(4, 4)

        sym_df = pd.DataFrame(df_mat)
        asym_df = sym_df.iloc[:3, :4]

        # Symmetric test_df, symmetric bg_df
        (is_test_df_sym1, is_bg_df_sym1) = sip.check_symmetry(sym_df, sym_df)
        self.assertTrue(is_test_df_sym1)
        self.assertTrue(is_bg_df_sym1)

        # Assymmetric test_df, symmetric bg_df
        (is_test_df_sym2, is_bg_df_sym2) = sip.check_symmetry(asym_df, sym_df)
        self.assertFalse(is_test_df_sym2)
        self.assertTrue(is_bg_df_sym2)

        # Assymetric bg should raise error
        with self.assertRaises(AssertionError) as e:
            sip.check_symmetry(sym_df, asym_df)
        self.assertIn("bg_df must be symmetric!", str(e.exception))

    def test_create_aggregated_fields_in_GCTs(self):

        # Make test_gct
        test_rids = ["M", "L", "P"]
        test_cids = ["Z", "X", "Y"]
        test_col_df = pd.DataFrame({"a": [1, 5, 6], "b": ["v", "l", "p"]})
        test_col_df.index = test_cids
        test_row_df = pd.DataFrame({"D": ["bee", "bird", "dog"],
                                    "C": ["bee", "me", "vee"]})
        test_row_df.index = test_rids
        test_gct = GCToo.GCToo(
            data_df=pd.DataFrame(np.nan, index=test_rids, columns=test_cids),
            row_metadata_df=test_row_df,
            col_metadata_df=test_col_df)

        # Make bg_gct
        bg_ids = ["u", "w", "v"]
        bg_meta_df = pd.DataFrame(index=bg_ids)
        bg_gct = GCToo.GCToo(
            data_df=pd.DataFrame(np.nan, index=bg_ids, columns=bg_ids),
            row_metadata_df=bg_meta_df,
            col_metadata_df=bg_meta_df.copy(deep=True))

        # Make expected results
        e_test_col_df = test_col_df.copy(deep=True)
        e_test_col_df2 = test_col_df.copy(deep=True)
        e_test_col_df["query_out"] = ["v|1", "l|5", "p|6"]
        e_test_col_df2["query_out"] = e_test_col_df2.index

        e_test_row_df = test_row_df.copy(deep=True)
        e_test_row_df["target_out"] = ["bee", "me", "vee"]

        e_bg_meta_df = bg_meta_df.copy(deep=True)
        e_bg_meta_df["target_out"] = ["u", "w", "v"]

        # Happy path
        out_test_gct, out_bg_gct = sip.create_aggregated_fields_in_GCTs(
            test_gct, bg_gct, ["b", "a"], ["C"], [], "query_out",
            "target_out", "|")

        pd.util.testing.assert_frame_equal(out_test_gct.col_metadata_df, e_test_col_df)
        pd.util.testing.assert_frame_equal(out_test_gct.row_metadata_df, e_test_row_df)
        pd.util.testing.assert_frame_equal(out_bg_gct.row_metadata_df, e_bg_meta_df)
        pd.util.testing.assert_frame_equal(out_bg_gct.col_metadata_df, e_bg_meta_df)

        # fields_to_aggregate_in_test_gct_queries is empty
        out_test_gct2, out_bg_gct2 = sip.create_aggregated_fields_in_GCTs(
            test_gct, bg_gct, [], ["C"], [], "query_out", "target_out", "|")

        pd.util.testing.assert_frame_equal(out_test_gct2.col_metadata_df, e_test_col_df2)
        pd.util.testing.assert_frame_equal(out_test_gct2.row_metadata_df, e_test_row_df)

    def test_aggregate_fields(self):
        df = pd.DataFrame({"a": ["a", "b", "c"],
                           "b": ["y", "l", "z"],
                           "c": [1, 6, 7]})
        out_col = ["a:1", "b:6", "c:7"]

        # Happy path
        out_df = sip.aggregate_fields(df, ["a", "c"], ":", "new_col")
        logger.debug("out_df:\n{}".format(out_df))

        df["new_col"] = out_col
        pd.util.testing.assert_frame_equal(out_df, df)

        # Metadata field supplied that's not actually present
        with self.assertRaises(AssertionError) as e:
            sip.aggregate_fields(df, ["d"], "blah", "blah")
        self.assertIn("d is not present", str(e.exception))

    def test_aggregate_metadata(self):
        df = pd.DataFrame({"pert_time": [24, 24, 24, 6, 6, 6],
                           "pert_id": ["A", "A", "A", "B", "B", "C"],
                           "pert_name": ["a", "A", "aa", "bee", "be", "B"],
                           "AGG": ["Y", "Y", "Y", "X", "X", "X"]})
        e_df = pd.DataFrame({"pert_time": ["6", "24"],
                             "pert_id": ["B|C", "A"],
                             "pert_name": ["B|be|bee", "A|a|aa"]})
        e_df.index = ["X", "Y"]

        out_df = sip.aggregate_metadata(df, "AGG", "|")
        logger.debug("out_df:\n{}".format(out_df))
        logger.debug("e_df:\n{}".format(e_df))
        pd.util.testing.assert_frame_equal(e_df, out_df, check_names=False)

        # Test a dataframe with just one sample
        e_df2 = pd.DataFrame([["A", "a", "24"]], index=["Y"],
                             columns=["pert_id", "pert_name", "pert_time"])
        out_df = sip.aggregate_metadata(df.iloc[[0], :], "AGG", "|")
        logger.debug("out_df:\n{}".format(out_df))
        pd.util.testing.assert_frame_equal(e_df2, out_df, check_names=False)

    def test_aggregate_one_series_uniquely(self):

        my_ser = pd.Series(["a", 3, 11])
        e_str = "3:11:a"

        out_str = sip.aggregate_one_series_uniquely(my_ser, sep=":")
        self.assertEqual(e_str, out_str)

    def test_extract_test_vals(self):

        # Symmetric GCT
        sym_test_data_df = pd.DataFrame(
            [[1.0, 0.5, 1.0, -0.4, 1.1, -0.6],
             [0.5, 1.0, 1.2, -0.8, -0.9, 0.4],
             [1.0, 1.2, 1.0, 0.1, 0.3, 1.3],
             [-0.4, -0.8, 0.1, 1.0, 0.5, -0.2],
             [1.1, -0.9, 0.3, 0.5, 1.0, 0.7],
             [-0.6, 0.4, 1.3, -0.2, 0.7, 1.0]])
        sym_test_meta_df = pd.DataFrame({
            "group": ["A", "B", "A", "B", "C", "C"],
            "id": [1, 2, 3, 4, 5, 6]})
        sym_test_gct = GCToo.GCToo(data_df=sym_test_data_df,
                                   row_metadata_df=sym_test_meta_df,
                                   col_metadata_df=sym_test_meta_df)

        # Expected values
        e_A_B_vals = [0.5, -0.4, 1.2, 0.1]
        e_A_C_vals = [1.1, 0.3, -0.6, 1.3]
        e_C_A_vals = [1.1, 0.3, -0.6, 1.3]
        e_A_A_vals = [1.0]

        A_B_vals = sip.extract_test_vals("A", "B", "group", "group", sym_test_gct, True)
        self.assertItemsEqual(e_A_B_vals, A_B_vals)

        A_C_vals = sip.extract_test_vals("A", "C", "group", "group", sym_test_gct, True)
        self.assertItemsEqual(e_A_C_vals, A_C_vals)

        C_A_vals = sip.extract_test_vals("C", "A", "group", "group", sym_test_gct, True)
        self.assertItemsEqual(e_C_A_vals, C_A_vals)

        A_A_vals = sip.extract_test_vals("A", "A", "group", "group", sym_test_gct, True)
        self.assertItemsEqual(e_A_A_vals, A_A_vals)

        # Verify that assert statement works
        with self.assertRaises(AssertionError) as e:
            sip.extract_test_vals("A", "D", "group", "group", sym_test_gct, True)
        self.assertIn("target D is not in the group metadata", str(e.exception))

        # Assymmetric GCT
        nonsym_test_row_meta_df = pd.DataFrame({
            "group": ["A", "B", "A", "B"],
            "id": [1, 2, 3, 4]})
        nonsym_test_col_meta_df = pd.DataFrame({
            "alt_group": ["F", "F", "E", "E"],
            "id": [1, 2, 3, 4]})
        nonsym_test_data_df = pd.DataFrame(
            [[1, 2, 3, 5],
             [7, 11, 13, 17],
             [19, 23, 29, 31],
             [-3, 5, 7, 11]])

        nonsym_test_gct = GCToo.GCToo(data_df=nonsym_test_data_df,
                                      row_metadata_df=nonsym_test_row_meta_df,
                                      col_metadata_df=nonsym_test_col_meta_df)

        # Expected values
        e_E_A_vals = [3, 5, 29, 31]
        e_F_B_vals = [7, 11, -3, 5]

        E_A_vals = sip.extract_test_vals("E", "A", "alt_group", "group", nonsym_test_gct, False)
        self.assertItemsEqual(e_E_A_vals, E_A_vals)

        F_B_vals = sip.extract_test_vals("F", "B", "alt_group", "group", nonsym_test_gct, False)
        self.assertItemsEqual(e_F_B_vals, F_B_vals)

    def test_extract_bg_vals_from_sym(self):

        bg_meta_df = pd.DataFrame({
            "group": ["A", "B", "A", "B", "C", "C"],
            "id": [1, 2, 3, 4, 5, 6]})
        bg_data_df = pd.DataFrame(
            [[1.0, 0.5, 1.0, -0.4, 1.1, -0.6],
             [0.5, 1.0, 1.2, -0.8, -0.9, 0.4],
             [1.0, 1.2, 1.0, 0.1, 0.3, 1.3],
             [-0.4, -0.8, 0.1, 1.0, 0.5, -0.2],
             [1.1, -0.9, 0.3, 0.5, 1.0, 0.7],
             [-0.6, 0.4, 1.3, -0.2, 0.7, 1.0]])
        bg_gct = GCToo.GCToo(data_df=bg_data_df,
                             row_metadata_df=bg_meta_df,
                             col_metadata_df=bg_meta_df)

        # Expected values
        e_A_vals = [0.5, 1.0, -0.4, 1.1, -0.6, 1.2, 0.1, 0.3, 1.3]
        e_B_vals = [0.5, 1.2, -0.8, -0.9, 0.4, -0.4, 0.1, 0.5, -0.2]
        e_C_vals = [1.1, -0.9, 0.3, 0.5, 0.7, -0.6, 0.4, 1.3, -0.2]

        A_vals = sip.extract_bg_vals_from_sym("A", "group", bg_gct)
        self.assertItemsEqual(e_A_vals, A_vals)

        B_vals = sip.extract_bg_vals_from_sym("B", "group", bg_gct)
        self.assertItemsEqual(e_B_vals, B_vals)

        C_vals = sip.extract_bg_vals_from_sym("C", "group", bg_gct)
        self.assertItemsEqual(e_C_vals, C_vals)

        # Verify that assert statement works
        with self.assertRaises(AssertionError) as e:
            sip.extract_bg_vals_from_sym("D", "group", bg_gct)
        self.assertIn("D is not in the group metadata", str(e.exception))

    def test_extract_bg_vals_from_non_sym(self):
        bg_row_meta_df = pd.DataFrame({
            "group": ["A", "B", "A", "B"],
            "id": [1, 2, 3, 4]})
        bg_col_meta_df = pd.DataFrame({
            "group": ["F", "F", "E", "E"],
            "id": [1, 2, 3, 4]})
        bg_data_df = pd.DataFrame(
            [[1, 2, 3, 5],
             [7, 11, 13, 17],
             [19, 23, 29, 31],
             [-3, 5, 7, 11]])
        bg_gct = GCToo.GCToo(data_df=bg_data_df,
                             row_metadata_df=bg_row_meta_df,
                             col_metadata_df=bg_col_meta_df)

        # Expected values
        e_A_vals = [1, 2, 3, 5, 19, 23, 29, 31]
        e_B_vals = [7, 11, 13, 17, -3, 5, 7, 11]

        A_vals = sip.extract_bg_vals_from_non_sym("A", "group", bg_gct)
        self.assertItemsEqual(e_A_vals, A_vals)

        B_vals = sip.extract_bg_vals_from_non_sym("B", "group", bg_gct)
        self.assertItemsEqual(e_B_vals, B_vals)

        # Verify that assert statement works
        with self.assertRaises(AssertionError) as e:
            sip.extract_bg_vals_from_non_sym("D", "group", bg_gct)
        self.assertIn("target D is not in the group metadata", str(e.exception))

    def test_percentile_score_single(self):
        test_vals = [7, 11, 13]
        bg_vals = [9, 11, -1, 19, 17, 7]

        out_score = sip.percentile_score_single(test_vals, bg_vals)
        self.assertAlmostEqual(out_score, 55.555, places=2)

    def test_compute_connectivities(self):

        # Create test_gct
        test_col_meta_df = pd.DataFrame({
            "pert": ["D", "D", "D", "E", "E", "E"],
            "cell": ["A375", "A375", "A375", "A375", "A375", "A375"],
            "agg": ["D:A375", "D:A375", "D:A375", "E:A375", "E:A375", "E:A375"],
            "other": ["M", "M", "N", "R", "P", "Q"],
            "other2": [3, 6, 4, 1, 1, 1.1]})

        test_row_meta_df = pd.DataFrame({
            "pert": ["A", "A", "B", "B"],
            "cell": ["A375", "A375", "A375", "A375"],
            "agg2": ["A:A375", "A:A375", "B:A375", "B:A375"],
            "weird": ["x", "y", "z", "z"]})

        test_data_df = pd.DataFrame(
            [[0.1, -0.3, -0.1, -0.4, 0.6, -0.7],
             [0.5, -0.7, -0.2, -1, 0.4, 0.2],
             [-0.2, 0.3, 0.7, 0.1, 0.4, -0.9],
             [0.1, 0.4, 0.2, 0.6, 0.4, -0.1]])
        test_gct = GCToo.GCToo(data_df=test_data_df,
                               row_metadata_df=test_row_meta_df,
                               col_metadata_df=test_col_meta_df)

        # Create bg_gct
        bg_meta_df = pd.DataFrame({
            "pert": ["A", "B", "A", "B", "C", "C"],
            "cell": ["A375", "A375", "A375", "A375", "A375", "A375"],
            "AGG": ["A:A375", "B:A375", "A:A375", "B:A375", "C:A375", "C:A375"],
            "ignore": ["j", "k", "l", "a", "b", "D"]})

        bg_data_df = pd.DataFrame(
            [[1.0, 0.5, 1.0, -0.4, 1.1, -0.6],
             [0.5, 1.0, 1.2, -0.8, -0.9, 0.4],
             [1.0, 1.2, 1.0, 0.1, 0.3, 1.3],
             [-0.4, -0.8, 0.1, 1.0, 0.5, -0.2],
             [1.1, -0.9, 0.3, 0.5, 1.0, 0.7],
             [-0.6, 0.4, 1.3, -0.2, 0.7, 1.0]])

        bg_gct = GCToo.GCToo(data_df=bg_data_df,
                             row_metadata_df=bg_meta_df,
                             col_metadata_df=bg_meta_df)

        # Create expected output
        A_bg = [0.5, 1.0, -0.4, 1.1, -0.6, 1.2, 0.1, 0.3, 1.3] # med = 0.4
        B_bg = [0.5, 1.2, -0.8, -0.9, 0.4, -0.4, 0.1, 0.5, -0.2] # med = 0.1
        (e_D_v_A, _) = stats.ks_2samp([0.1, -0.3, -0.1, 0.5, -0.7, -0.2], A_bg) # med = -1.5, so -
        (e_D_v_B, _) = stats.ks_2samp([-0.2, 0.3, 0.7, 0.1, 0.4, 0.2], B_bg) # med = 0.25, so +
        (e_E_v_A, _) = stats.ks_2samp([-0.4, 0.6, -0.7, -1, 0.4, 0.2], A_bg) # med = -0.1, so -
        (e_E_v_B, _) = stats.ks_2samp([0.1, 0.4, -0.9, 0.6, 0.4, -0.1], B_bg) # med = 0.25, so +

        e_conn_df = pd.DataFrame(
            [[e_D_v_A, e_E_v_A], [e_D_v_B, e_E_v_B]],
            index = ["A:A375", "B:A375"],
            columns = ["D:A375", "E:A375"])
        e_signed_conn_df = pd.DataFrame(
            [[-e_D_v_A, -e_E_v_A], [e_D_v_B, e_E_v_B]],
            index = ["A:A375", "B:A375"],
            columns = ["D:A375", "E:A375"])

        e_row_meta_df = pd.DataFrame({
            "pert": ["A", "B"],
            "cell": ["A375", "A375"]})
        e_row_meta_df.index = ["A:A375", "B:A375"]

        e_row_meta_df = pd.DataFrame({
            "pert": ["A", "B"],
            "cell": ["A375", "A375"],
            "weird": ["x:y", "z"]})
        e_row_meta_df.index = ["A:A375", "B:A375"]

        e_col_meta_df = pd.DataFrame({
            "pert": ["D", "E"],
            "cell": ["A375", "A375"],
            "other": ["M:N", "P:Q:R"],
            "other2": ["3.0:4.0:6.0", "1.0:1.1"]})
        e_col_meta_df.index = ["D:A375", "E:A375"]

        (conn_gct, signed_conn_gct) = sip.compute_connectivities(
            test_gct, bg_gct, "agg", "agg2", "AGG", "ks_test", False, ":")

        logger.debug("conn_gct.data_df:\n{}".format(conn_gct.data_df))
        logger.debug("e_conn_df:\n{}".format(e_conn_df))
        logger.debug("conn_gct.row_metadata_df:\n{}".format(conn_gct.row_metadata_df))
        logger.debug("conn_gct.col_metadata_df:\n{}".format(conn_gct.col_metadata_df))

        pd.util.testing.assert_frame_equal(conn_gct.data_df, e_conn_df)
        pd.util.testing.assert_frame_equal(signed_conn_gct.data_df, e_signed_conn_df)
        pd.util.testing.assert_frame_equal(conn_gct.row_metadata_df, e_row_meta_df, check_names=False)
        pd.util.testing.assert_frame_equal(conn_gct.col_metadata_df, e_col_meta_df, check_names=False)

        # Make sure connectivity metric is valid
        with self.assertRaises(Exception) as e:
            sip.compute_connectivities(test_gct, bg_gct, "agg",
                                       "agg2", "AGG", "wtcs",
                                       False, "|")
        self.assertIn("connectivity metric must be either ks_test or", str(e.exception))

        # Make sure we have agreement across test_gct and bg_gct
        with self.assertRaises(Exception) as e:
            sip.compute_connectivities(test_gct, bg_gct, "agg", "pert",
                                       "ignore", "wtcs", False, "|")
        self.assertIn("There are no targets ", str(e.exception))


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()