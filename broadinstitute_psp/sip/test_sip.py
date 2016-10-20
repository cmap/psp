import unittest
import logging
import os
import numpy as np
import pandas as pd
import scipy.stats as stats

import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_cmap.io.GCToo.parse_gctoo as pg
import sip

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

FUNCTIONAL_TESTS_DIR = "sip/functional_tests"


class TestSip(unittest.TestCase):

    def test_main(self):
        test_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_sip_in_test.gct")
        bg_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_sip_in_bg.gct")
        out_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_sip_main_out.gct")

        args_string = "-t {} -b {} -o {} -tfq {} -tft {} -bf {}".format(
            test_gct_path, bg_gct_path, out_path, "pert_iname", "pert_iname", "pert_iname")
        args = sip.build_parser().parse_args(args_string.split())

        # Run main method
        sip.main(args)

        # Compare the output of main with the expected output
        e_out_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_sip_expected_conn.gct")
        e_out_gct = pg.parse(e_out_path)
        out_gct = pg.parse(out_path)

        self.assertTrue(np.allclose(e_out_gct.data_df.values, out_gct.data_df.values), (
            "\ne_out_gct.data_df:\n{}\nout_gct.data_df:\n{}".format(
                e_out_gct.data_df, out_gct.data_df)))
        self.assertTrue(e_out_gct.row_metadata_df.equals(out_gct.row_metadata_df), (
            "\ne_out_gct.row_metadata_df:\n{}\nout_gct.row_metadata_df:\n{}".format(
                e_out_gct.row_metadata_df, out_gct.row_metadata_df)))
        self.assertTrue(e_out_gct.col_metadata_df.equals(out_gct.col_metadata_df), (
            "\ne_out_gct.col_metadata_df:\n{}\nout_gct.col_metadata_df:\n{}".format(
                e_out_gct.col_metadata_df, out_gct.col_metadata_df)))

        # Remove the created file
        os.remove(out_path)

    def test_prepare_multi_index_dfs(self):

        test_df_index = pd.MultiIndex.from_arrays(
            [["A", "B", "C", "C", "B", "A"],
             ["A375", "A375", "A375", "A375", "A375", "A375"],
             [10, 10, 10, 10, 10, 10]], names=["pert_iname", "cell", "dose"])
        test_df_columns = pd.MultiIndex.from_arrays(
            [["F", "F", "F", "E", "E", "E"], [1, 2, 3, 4, 5, 6]], names=["pert_uname", "dose"])
        test_df = pd.DataFrame(
            [[1, 3, 5, 7, 11.0, 13.0],
             [17, 19, 23, 29, 31, 37],
             [-1.0, -1.1, -1.3, -1.5, 1.7, -1.9],
             [11, 9, 7, 5, 3, -1],
             [1.0, -0.9, 0.8, -0.7, 0.6, -0.5],
             [-77, 66.0, 55.1, 44.3, 51.2, 1.0]],
            index=test_df_index, columns=test_df_columns)

        bg_df_index = pd.MultiIndex.from_arrays(
            [["E123", "F21", "C666", "D1"], [1, 2, 3, 4]], names=["cell", "thing"])
        bg_df = pd.DataFrame(
            [[1, 2, 3, 5],
             [7, 11, 13, 17],
             [19, 23, 29, 31],
             [-3, 5, 7, 11]],
            index=bg_df_index, columns=bg_df_index)

        e_test_df_index = pd.MultiIndex.from_arrays(
            [["A", "A", "B", "B", "C", "C"],
             [10, 10, 10, 10, 10, 10],
             ["A:10", "A:10", "B:10", "B:10", "C:10", "C:10"]], names=["pert_iname", "dose", "target_field"])
        e_test_df_columns = pd.MultiIndex.from_arrays(
            [["E", "E", "E", "F", "F", "F"], ["E", "E", "E", "F", "F", "F"]], names=["pert_uname", "query_field"])
        e_test_df = pd.DataFrame(
            [[7, 11.0, 13.0, 1, 3, 5],
             [44.3, 51.2, 1.0, -77, 66.0, 55.1],
             [29, 31, 37, 17, 19, 23],
             [-0.7, 0.6, -0.5, 1.0, -0.9, 0.8],
             [-1.5, 1.7, -1.9, -1.0, -1.1, -1.3],
             [5, 3, -1, 11, 9, 7]],
            index=e_test_df_index, columns=e_test_df_columns)
        e_bg_df_index = pd.MultiIndex.from_arrays(
            [["C666", "D1", "E123", "F21"], ["C666", "D1", "E123", "F21"]], names=["cell", "target_field"])
        e_bg_df = pd.DataFrame(
            [[29, 31, 19, 23],
             [7, 11, -3, 5],
             [3, 5, 1, 2],
             [13, 17, 7, 11]],
            index=e_bg_df_index, columns=e_bg_df_index)

        (out_test_df, out_bg_df) = sip.prepare_multi_index_dfs(
            test_df, bg_df, ["pert_uname"], ["pert_iname", "dose"], ["cell"], "query_field", "target_field", ":")

        pd.util.testing.assert_frame_equal(out_test_df, e_test_df)
        pd.util.testing.assert_frame_equal(out_bg_df, e_bg_df)

    def test_check_symmetry(self):
        df_mat = np.array(
            [[29, 31, 19, 23],
             [7, 11, -3, 5],
             [3, 5, 1, 2],
             [13, 17, 7, 11]])

        # Asymmetric test_df
        test_df_index = pd.MultiIndex.from_arrays(
            [["A", "B", "C", "D"]], names=["sample_id", "plate_id"])
        asym_test_df_columns = pd.MultiIndex.from_arrays(
            [["A", "B", "C", "D"], ["X1", "X2", "X3", "X4"]],
            names=["sample_id", "plate_id"])
        asym_test_df = pd.DataFrame(df_mat, index=test_df_index, columns=asym_test_df_columns)

        # Symmetric test_df
        sym_test_df = pd.DataFrame(df_mat, index=test_df_index, columns=test_df_index)

        # Asymmetric bg_df
        sym_bg_df_index = pd.MultiIndex.from_arrays(
            [["A", "B", "C", "D"]], names=["cell"])
        asym_bg_df_columns = pd.MultiIndex.from_arrays(
            [["D", "E", "F", "G"]], names=["cell"])
        asym_bg_df = pd.DataFrame(df_mat, index=sym_bg_df_index, columns=asym_bg_df_columns)

        # Symmetric bg_df
        sym_bg_df = pd.DataFrame(df_mat, index=sym_bg_df_index, columns=sym_bg_df_index)

        # Symmetric test_df, symmetric bg_df
        (is_test_df_sym1, is_bg_df_sym1) = sip.check_symmetry(sym_test_df, sym_bg_df)
        self.assertTrue(is_test_df_sym1)
        self.assertTrue(is_bg_df_sym1)

        # Assymmetric test_df, symmetric bg_df
        (is_test_df_sym2, is_bg_df_sym2) = sip.check_symmetry(asym_test_df, sym_bg_df)
        self.assertFalse(is_test_df_sym2)
        self.assertTrue(is_bg_df_sym2)

        # Assymetric bg should raise error
        with self.assertRaises(AssertionError) as e:
            sip.check_symmetry(sym_test_df, asym_bg_df)
        self.assertIn("bg_mi_df must be symmetric!", str(e.exception))

    def test_extract_test_vals(self):
        # Symmetric
        sym_test_df_index = pd.MultiIndex.from_arrays(
            [["A", "B", "A", "B", "C", "C"], [1, 2, 3, 4, 5, 6]], names=["group", "id"])
        sym_test_df = pd.DataFrame(
            [[1.0, 0.5, 1.0, -0.4, 1.1, -0.6],
             [0.5, 1.0, 1.2, -0.8, -0.9, 0.4],
             [1.0, 1.2, 1.0, 0.1, 0.3, 1.3],
             [-0.4, -0.8, 0.1, 1.0, 0.5, -0.2],
             [1.1, -0.9, 0.3, 0.5, 1.0, 0.7],
             [-0.6, 0.4, 1.3, -0.2, 0.7, 1.0]],
            index=sym_test_df_index, columns=sym_test_df_index)

        # Expected values
        e_A_B_vals = [0.5, -0.4, 1.2, 0.1]
        e_A_C_vals = [1.1, 0.3, -0.6, 1.3]
        e_C_A_vals = [1.1, 0.3, -0.6, 1.3]
        e_A_A_vals = [1.0]

        A_B_vals = sip.extract_test_vals("A", "B", "group", "group", sym_test_df, True)
        self.assertItemsEqual(e_A_B_vals, A_B_vals)

        A_C_vals = sip.extract_test_vals("A", "C", "group", "group", sym_test_df, True)
        self.assertItemsEqual(e_A_C_vals, A_C_vals)

        C_A_vals = sip.extract_test_vals("C", "A", "group", "group", sym_test_df, True)
        self.assertItemsEqual(e_C_A_vals, C_A_vals)

        A_A_vals = sip.extract_test_vals("A", "A", "group", "group", sym_test_df, True)
        self.assertItemsEqual(e_A_A_vals, A_A_vals)

        # Verify that assert statement works
        with self.assertRaises(AssertionError) as e:
            sip.extract_test_vals("A", "D", "group", "group", sym_test_df, True)
        self.assertIn("target D is not in the group level", str(e.exception))

        # Assymmetric
        nonsym_test_df_index = pd.MultiIndex.from_arrays(
            [["A", "B", "A", "B"], [1, 2, 3, 4]], names=["group", "id"])
        nonsym_test_df_columns = pd.MultiIndex.from_arrays(
            [["F", "F", "E", "E"], [1, 2, 3, 4]], names=["alt_group", "id"])
        nonsym_test_df = pd.DataFrame(
            [[1, 2, 3, 5],
             [7, 11, 13, 17],
             [19, 23, 29, 31],
             [-3, 5, 7, 11]],
            index=nonsym_test_df_index, columns=nonsym_test_df_columns)

        # Expected values
        e_E_A_vals = [3, 5, 29, 31]
        e_F_B_vals = [7, 11, -3, 5]

        E_A_vals = sip.extract_test_vals("E", "A", "alt_group", "group", nonsym_test_df, False)
        self.assertItemsEqual(e_E_A_vals, E_A_vals)

        F_B_vals = sip.extract_test_vals("F", "B", "alt_group", "group", nonsym_test_df, False)
        self.assertItemsEqual(e_F_B_vals, F_B_vals)

    def test_extract_bg_vals_from_sym(self):
        bg_df_index = pd.MultiIndex.from_arrays(
            [["A", "B", "A", "B", "C", "C"], [1, 2, 3, 4, 5, 6]], names=["group", "id"])
        bg_df = pd.DataFrame(
            [[1.0, 0.5, 1.0, -0.4, 1.1, -0.6],
             [0.5, 1.0, 1.2, -0.8, -0.9, 0.4],
             [1.0, 1.2, 1.0, 0.1, 0.3, 1.3],
             [-0.4, -0.8, 0.1, 1.0, 0.5, -0.2],
             [1.1, -0.9, 0.3, 0.5, 1.0, 0.7],
             [-0.6, 0.4, 1.3, -0.2, 0.7, 1.0]],
            index=bg_df_index, columns=bg_df_index)

        # Expected values
        e_A_vals = [0.5, 1.0, -0.4, 1.1, -0.6, 1.2, 0.1, 0.3, 1.3]
        e_B_vals = [0.5, 1.2, -0.8, -0.9, 0.4, -0.4, 0.1, 0.5, -0.2]
        e_C_vals = [1.1, -0.9, 0.3, 0.5, 0.7, -0.6, 0.4, 1.3, -0.2]

        A_vals = sip.extract_bg_vals_from_sym("A", "group", bg_df)
        self.assertItemsEqual(e_A_vals, A_vals)

        B_vals = sip.extract_bg_vals_from_sym("B", "group", bg_df)
        self.assertItemsEqual(e_B_vals, B_vals)

        C_vals = sip.extract_bg_vals_from_sym("C", "group", bg_df)
        self.assertItemsEqual(e_C_vals, C_vals)

        # Verify that assert statement works
        with self.assertRaises(AssertionError) as e:
            sip.extract_bg_vals_from_sym("D", "group", bg_df)
        self.assertIn("D is not in the group level", str(e.exception))

    def test_extract_bg_vals_from_non_sym(self):
        bg_df_index = pd.MultiIndex.from_arrays(
            [["A", "B", "A", "B"], [1, 2, 3, 4]], names=["group", "id"])
        bg_df_columns = pd.MultiIndex.from_arrays(
            [["F", "F", "E", "E"], [1, 2, 3, 4]], names=["group", "id"])
        bg_df = pd.DataFrame(
            [[1, 2, 3, 5],
             [7, 11, 13, 17],
             [19, 23, 29, 31],
             [-3, 5, 7, 11]],
            index=bg_df_index, columns=bg_df_columns)

        # Expected values
        e_A_vals = [1, 2, 3, 5, 19, 23, 29, 31]
        e_B_vals = [7, 11, 13, 17, -3, 5, 7, 11]

        A_vals = sip.extract_bg_vals_from_non_sym("A", "group", bg_df)
        self.assertItemsEqual(e_A_vals, A_vals)

        B_vals = sip.extract_bg_vals_from_non_sym("B", "group", bg_df)
        self.assertItemsEqual(e_B_vals, B_vals)

        # Verify that assert statement works
        with self.assertRaises(AssertionError) as e:
            sip.extract_bg_vals_from_non_sym("D", "group", bg_df)
        self.assertIn("target D is not in the group level", str(e.exception))

    def test_percentile_score_single(self):
        test_vals = [7, 11, 13]
        bg_vals = [9, 11, -1, 19, 17, 7]

        out_score = sip.percentile_score_single(test_vals, bg_vals)
        self.assertAlmostEqual(out_score, 55.555, places=2)

    def test_compute_connectivities(self):
        # External query against build
        test_df_index = pd.MultiIndex.from_arrays(
            [["A", "A", "B", "B"], ["A375", "A375", "A375", "A375"],
             ["A:A375", "A:A375", "B:A375", "B:A375"]], names=["pert", "cell", "aggregated"])
        test_df_columns = pd.MultiIndex.from_arrays(
            [["D", "D", "D", "E", "E", "E"], ["A375", "A375", "A375", "A375", "A375", "A375"],
             ["D:A375", "D:A375", "D:A375", "E:A375", "E:A375", "E:A375"]],
            names=["pert_iname", "cell", "aggregated2"])
        test_df = pd.DataFrame(
            [[0.1, -0.3, -0.1, -0.4, 0.6, -0.7],
             [0.5, -0.7, -0.2, -1, 0.4, 0.2],
             [-0.2, 0.3, 0.7, 0.1, 0.4, -0.9],
             [0.1, 0.4, 0.2, 0.6, 0.4, -0.1]],
            index=test_df_index, columns=test_df_columns)

        bg_df_index = pd.MultiIndex.from_arrays(
            [["A", "B", "A", "B", "C", "C"], ["A375", "A375", "A375", "A375", "A375", "A375"],
             ["A:A375", "B:A375", "A:A375", "B:A375", "C:A375", "C:A375"]],
            names=["pert", "cell", "bg_aggregated"])
        bg_df = pd.DataFrame(
            [[1.0, 0.5, 1.0, -0.4, 1.1, -0.6],
             [0.5, 1.0, 1.2, -0.8, -0.9, 0.4],
             [1.0, 1.2, 1.0, 0.1, 0.3, 1.3],
             [-0.4, -0.8, 0.1, 1.0, 0.5, -0.2],
             [1.1, -0.9, 0.3, 0.5, 1.0, 0.7],
             [-0.6, 0.4, 1.3, -0.2, 0.7, 1.0]],
            index=bg_df_index, columns=bg_df_index)

        A_bg = [0.5, 1.0, -0.4, 1.1, -0.6, 1.2, 0.1, 0.3, 1.3] # med = 0.4
        B_bg = [0.5, 1.2, -0.8, -0.9, 0.4, -0.4, 0.1, 0.5, -0.2] # med = 0.1
        (e_D_v_A, _) = stats.ks_2samp([0.1, -0.3, -0.1, 0.5, -0.7, -0.2], A_bg) # med = -1.5, so -
        (e_D_v_B, _) = stats.ks_2samp([-0.2, 0.3, 0.7, 0.1, 0.4, 0.2], B_bg) # med = 0.25, so +
        (e_E_v_A, _) = stats.ks_2samp([-0.4, 0.6, -0.7, -1, 0.4, 0.2], A_bg) # med = -0.1, so -
        (e_E_v_B, _) = stats.ks_2samp([0.1, 0.4, -0.9, 0.6, 0.4, -0.1], B_bg) # med = 0.25, so +

        e_conn_df_index = pd.MultiIndex.from_arrays(
            [["A", "B"], ["A375", "A375"], ["A:A375", "B:A375"]],
            names=["pert", "cell", "aggregated"])
        e_conn_df_columns = pd.MultiIndex.from_arrays(
            [["D", "E"], ["A375", "A375"], ["D:A375", "E:A375"]],
            names=["pert_iname", "cell", "aggregated2"])
        e_conn_df = pd.DataFrame(
            [[e_D_v_A, e_E_v_A], [e_D_v_B, e_E_v_B]], index=e_conn_df_index, columns=e_conn_df_columns)
        e_signed_conn_df = pd.DataFrame(
            [[-e_D_v_A, -e_E_v_A], [e_D_v_B, e_E_v_B]], index=e_conn_df_index, columns=e_conn_df_columns)

        (conn_df, signed_conn_df) = sip.compute_connectivities(
            test_df, bg_df, "aggregated2", "aggregated", "bg_aggregated", "ks_test", False)

        pd.util.testing.assert_frame_equal(conn_df, e_conn_df, (
            "\nconn_df:\n{}\ne_conn_df:\n{}").format(conn_df, e_conn_df))
        pd.util.testing.assert_frame_equal(signed_conn_df, e_signed_conn_df, (
            "\nsigned_conn_df:\n{}\ne_signed_conn_df:\n{}").format(
            signed_conn_df, e_signed_conn_df))

        # Check that assertion works
        with self.assertRaises(Exception) as e:
            sip.compute_connectivities(test_df, bg_df, "aggregated2", "aggregated", "bg_aggregated", "wtcs", False)
        self.assertIn("connectivity metric must be either ks_test or", str(e.exception))

    def test_add_aggregated_level_to_multi_index(self):
        # Create group_id column
        mi = pd.MultiIndex.from_arrays(
            [["a", "b", "c"], ["10", "10", "1"], ["24", "24", "24"]],
             names=["cell", "dose", "time"])
        e_out_mi = pd.MultiIndex.from_arrays(
            [["a", "b", "c"], ["10", "10", "1"], ["24", "24", "24"], ["a:10", "b:10", "c:1"]],
             names=["cell", "dose", "time", "aggregated"])
        e_out_subset_mi = pd.MultiIndex.from_arrays(
            [["a", "b", "c"], ["10", "10", "1"], ["a:10", "b:10", "c:1"]],
             names=["cell", "dose", "aggregated"])

        (out_mi, out_subset_mi) = sip.add_aggregated_level_to_multi_index(
            mi, ["cell", "dose"], "aggregated", sep=":")

        self.assertTrue(out_mi.equals(e_out_mi), (
            "\nout_mi:\n{}\ne_out_mi:\n{}".format(out_mi, e_out_mi)))
        self.assertTrue(out_subset_mi.equals(e_out_subset_mi), (
            "\nout_subset_mi:\n{}\ne_out_subset_mi:\n{}".format(out_subset_mi, e_out_subset_mi)))

        e_out_mi2 = pd.MultiIndex.from_arrays(
            [["a", "b", "c"], ["10", "10", "1"], ["24", "24", "24"], ["a", "b", "c"]],
             names=["cell", "dose", "time", "aggregated"])
        e_out_subset_mi2 = pd.MultiIndex.from_arrays(
            [["a", "b", "c"], ["a", "b", "c"]],
             names=["cell", "aggregated"])

        (out_mi2, out_subset_mi2) = sip.add_aggregated_level_to_multi_index(
            mi, ["cell"], "aggregated", sep=":")

        pd.util.testing.assert_index_equal(out_mi2, e_out_mi2)
        pd.util.testing.assert_index_equal(out_subset_mi2, e_out_subset_mi2)


if __name__ == "__main__":
    setup_logger.setup(verbose=False)
    unittest.main()