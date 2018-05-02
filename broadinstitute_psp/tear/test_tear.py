import unittest
import logging
import os
import numpy as np
import pandas as pd

import cmapPy.pandasGEXpress as GCToo
import cmapPy.pandasGEXpress.parse as parse
import broadinstitute_psp.utils.setup_logger as setup_logger
import tear

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Use functional tests assets from the tear directory
FUNCTIONAL_TESTS_DIR = "tear/functional_tests"


class TestTear(unittest.TestCase):
    def test_main(self):
        in_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_tear_main.gct")
        out_name = os.path.join(FUNCTIONAL_TESTS_DIR, "test_tear_out.gct")

        args_string = ("-i {} -o {} -dm -p {}").format(
            in_gct_path, out_name, "psp_production.cfg")
        args = tear.build_parser().parse_args(args_string.split())
        tear.main(args)

        # Read in result
        out_gct = parse.parse(out_name)

        e_values = np.array(
            [[0., 4.07, -1.48, -10.71, 0.],
            [4.43, -3.26, -0.23, 0., 1.48],
            [0., 2.49, 2.50, -1.48, -0.86]])
        self.assertTrue(np.allclose(e_values, out_gct.data_df, atol=1e-2))

        # Clean up
        os.remove(out_name)

    def test_median_normalize(self):
        data = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9], [10, 11, 12]],
                            index=["a", "b", "c", "d"],
                            columns=["e", "f", "g"],
                            dtype=float)
        row_meta = pd.DataFrame([["1", "rm2"],["1", "rm4"],["1", "rm6"],["1", "rm8"]],
                                index=["a", "b", "c", "d"],
                                columns=["row_subset", "row_field2"])
        col_meta = pd.DataFrame([["cm1", "1"],["cm3", "2"],["cm5", "2"]],
                                index=["e", "f", "g"],
                                columns=["col_field1", "col_subset"])
        in_gct = GCToo.GCToo.GCToo(data_df=data, row_metadata_df=row_meta, col_metadata_df=col_meta)
        config_metadata = {
            "row_subset_field": "row_subset",
            "col_subset_field": "col_subset",
            "subset_normalize_prov_code_entry": "GMN",
            "row_normalize_prov_code_entry": "RMN"}

        # Subset normalize
        ignore_subset_norm = False
        e_data = pd.DataFrame([[0, -0.5, 0.5], [0, -0.5, 0.5], [0, -0.5, 0.5], [0, -0.5, 0.5]],
                            index=["a", "b", "c", "d"],
                            columns=["e", "f", "g"])
        e_prov_code = ["A", "B", "GMN"]

        (out_gct, out_prov_code) = tear.median_normalize(
            in_gct, False, ignore_subset_norm, config_metadata, ["A", "B"])

        self.assertTrue(np.allclose(out_gct.data_df, e_data, atol=1e-2))
        self.assertEqual(out_prov_code, e_prov_code)

        # Ordinary normalization
        ignore_subset_norm2 = True
        e_data2 = pd.DataFrame([[-1, 0, 1], [-1, 0, 1], [-1, 0, 1], [-1, 0, 1]],
                               index=["a", "b", "c", "d"],
                               columns=["e", "f", "g"],
                               dtype=float)
        row_meta2 = pd.DataFrame([["1", "rm2"],["1", "rm4"],["1", "rm6"],["1", "rm8"]],
                                index=["a", "b", "c", "d"],
                                columns=["row_subset", "row_field2"])
        col_meta2 = pd.DataFrame([["cm1", "1"],["cm3", "1"],["cm5", "1"]],
                                index=["e", "f", "g"],
                                columns=["col_field1", "col_subset"])
        in_gct2 = GCToo.GCToo.GCToo(data_df=data, row_metadata_df=row_meta2, col_metadata_df=col_meta2)
        e_prov_code2 = ["A", "B", "RMN"]

        (out_gct2, out_prov_code2) = tear.median_normalize(
            in_gct2, False, ignore_subset_norm2, config_metadata, ["A", "B"])

        self.assertTrue(np.allclose(out_gct2.data_df, e_data2, atol=1e-2))
        self.assertEqual(out_prov_code2, e_prov_code2)

    def test_check_for_subsets(self):
        row_meta = pd.DataFrame([["rm1", "rm2"],["rm3", "rm4"],
                                 ["rm5", "rm6"],["rm7", "rm8"]],
                                index=["a", "b", "c", "d"],
                                columns=["row_field1", "row_field2"])
        col_meta = pd.DataFrame([["cm1", "cm2"],["cm3", "cm4"],["cm5", "cm6"]],
                                index=["e", "f", "g"],
                                columns=["col_field1", "col_field2"])
        row_field = "row_field1"
        col_field = "col_field2"
        subsets_exist = tear.check_for_subsets(row_meta, col_meta, row_field, col_field)

        self.assertTrue(subsets_exist)

    def test_row_median_normalize(self):
        df = pd.DataFrame(np.array([[10, -3, 1.2, 0.6],
                                    [0.45, 0.2, 0, 0.2],
                                    [4.5, np.nan, 0.3, 0.4]], dtype=float))
        e_df = pd.DataFrame(np.array([[9.1, -3.9, 0.3, -0.3],
                                    [0.25, 0, -0.2, 0],
                                    [4.1, np.nan, -0.1, 0]], dtype=float))
        out_df = tear.row_median_normalize(df, False)
        self.assertTrue(np.allclose(out_df, e_df, atol=1e-2, equal_nan=True),
                        ("\nExpected out_df:\n{} " +
                         "\nActual out_df:\n{}").format(e_df, out_df))

        # Divide by MAD this time
        out_df_mad = tear.row_median_normalize(df, True)
        e_df_mad = pd.DataFrame(np.array([
            [6.42, -2.75, 0.21, -0.21],
            [3.71, 0.00, -2.97, 0.00],
            [60.79, np.nan, -1.48, 0.00]], dtype=float))
        self.assertTrue(np.allclose(out_df_mad, e_df_mad, atol=1e-2, equal_nan=True),
                        ("\nExpected out_df:\n{} " +
                         "\nActual out_df:\n{}").format(e_df_mad, out_df_mad))

    def test_subset_normalize(self):
        ROW_SUBSET_FIELD = "pr_probe_normalization_group"
        COL_SUBSET_FIELD = "det_normalization_group_vector"
        row_df = pd.DataFrame(np.array([["8350", "1"], ["8350", "1"],
                                        ["8350", "2"], ["8350", "2"]]),
                              index=["r1", "r2", "r3", "r4"],
                              columns=["pr_gene_id", "pr_probe_normalization_group"])
        col_df = pd.DataFrame(np.array([["G-0022", "1,1"], ["G-0022", "1,1"], ["G-0022", "1,2"],
                                        ["G-0022", "2,2"], ["G-0022", "2,2"]]),
                              index=["c1", "c2", "c3", "c4", "c5"],
                              columns=["det_plate", "det_normalization_group_vector"])
        data_df = pd.DataFrame([[7, 8, 3, 8, 9],
                                [9, 7, 4, 9, 2],
                                [8, 6, 7, 8, 2],
                                [4, 8, 5, 5, 7]],
                               index=["r1", "r2", "r3", "r4"],
                               columns=["c1", "c2", "c3", "c4", "c5"],
                               dtype=float)
        in_gct = GCToo.GCToo.GCToo(data_df=data_df, row_metadata_df=row_df, col_metadata_df=col_df)
        e_df = pd.DataFrame(np.array([[0, 1, -4, -0.5, 0.5],
                                      [2, 0, -3, 3.5, -3.5],
                                      [1, -1, 0, 1, -5],
                                      [-2, 2, 0, 0, 2]], dtype=float))
        out_df = tear.subset_normalize(in_gct, False, ROW_SUBSET_FIELD, COL_SUBSET_FIELD)
        self.assertTrue(np.array_equal(out_df, e_df),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_df, out_df))

    def test_make_norm_ndarray(self):
        ROW_SUBSET_FIELD = "pr_probe_normalization_group"
        COL_SUBSET_FIELD = "det_normalization_group_vector"
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
        norm_ndarray = tear.make_norm_ndarray(row_df, col_df, ROW_SUBSET_FIELD, COL_SUBSET_FIELD)
        self.assertTrue(np.array_equal(norm_ndarray, e_norm_ndarray),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_norm_ndarray, norm_ndarray))

    def test_iterate_over_norm_ndarray_and_normalize(self):
        data_df = pd.DataFrame(np.array([[7, 8, 3, 8, 9],
                                         [9, 7, 4, 9, 2],
                                         [8, 6, 7, 8, 2],
                                         [4, 8, 5, 5, 7]]),
                                        index = ["a", "b", "c", "d"],
                                        dtype="float")
        norm_ndarray = np.array([[1, 1, 1, 2, 2],
                                 [1, 1, 1, 2, 2],
                                 [1, 1, 2, 2, 2],
                                 [1, 1, 2, 2, 2]])
        e_df = pd.DataFrame(np.array([[0, 1, -4, -0.5, 0.5],
                                      [2, 0, -3, 3.5, -3.5],
                                      [1, -1, 0, 1, -5],
                                      [-2, 2, 0, 0, 2]], dtype="float"))
        out_df = tear.iterate_over_norm_ndarray_and_normalize(data_df, norm_ndarray, False)
        self.assertTrue(np.array_equal(out_df, e_df),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_df, out_df))

        # Divide by MAD this time
        out_df_mad = tear.iterate_over_norm_ndarray_and_normalize(data_df, norm_ndarray, True)
        e_df_mad = pd.DataFrame(np.array([[0, 1.48, -5.93, -1.48, 1.48],
                                          [1.48, 0, -2.22, 1.48, -1.48],
                                          [1.48, -1.48, 0, 1.48, -7.41],
                                          [-1.48, 1.48, 0, 0, 2965159.37]], dtype="float"))
        self.assertTrue(np.allclose(out_df_mad, e_df_mad, atol=1e-2),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_df_mad, out_df_mad))

        # Slightly different norm_ndarray
        norm_ndarray = np.array([[1, 1, 1, 2, 2],
                                 [1, 1, 1, 2, 2],
                                 [1, 1, 2, 2, 3],
                                 [1, 1, 2, 2, 3]])
        e_df = pd.DataFrame(np.array([[0, 1, -4, -0.5, 0.5],
                                      [2, 0, -3, 3.5, -3.5],
                                      [1, -1, -0.5, 0.5, 0],
                                      [-2, 2, 0, 0, 0]], dtype="float"))
        out_df = tear.iterate_over_norm_ndarray_and_normalize(data_df, norm_ndarray, False)
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
        out_df = tear.iterate_over_norm_ndarray_and_normalize(data_df, norm_ndarray, False)
        self.assertTrue(np.array_equal(out_df, e_df),
                        ("\nExpected out:\n{} " +
                         "\nActual out:\n{}").format(e_df, out_df))


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
