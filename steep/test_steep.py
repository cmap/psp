import unittest
import logging
import utils.setup_logger as setup_logger
import os
import numpy as np
import pandas as pd
import scipy.stats as stats

import in_out.parse_gctoo as parse_gctoo
import in_out.GCToo as GCToo
import steep

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestSteep(unittest.TestCase):

    def test_compute_similarity_matrix(self):
        data_df = pd.DataFrame.from_dict(
            {"a":[3,5,1,2,4], "b":[2,9,np.nan,8,1], "c":[np.nan,3,5,2,8]})

        # Spearman
        spearman_out = steep.compute_similarity_matrix(data_df, method="spearman")

        # Confirm correlation that the correct items were used to compute corr
        a_v_b = pd.Series([3, 5, 2, 4]).corr(pd.Series([2, 9, 8, 1]), method="spearman")
        self.assertAlmostEqual(a_v_b, spearman_out.iloc[0, 1], msg=(
            "Expected value is {}. spearman_out.iloc[0, 1]: {}".format(
                a_v_b, spearman_out.iloc[0, 1])))
        a_v_c = pd.Series([5, 1, 2, 4]).corr(pd.Series([3, 5, 2, 8]), method="spearman")
        self.assertAlmostEqual(a_v_c, spearman_out.iloc[0, 2], msg=(
            "Expected value is {}. spearman_out.iloc[0, 2]: {}".format(
                a_v_c, spearman_out.iloc[0, 2])))
        b_v_c = pd.Series([9, 8, 1]).corr(pd.Series([3, 2, 8]), method="spearman")
        self.assertAlmostEqual(b_v_c, spearman_out.iloc[1, 2], msg=(
            "Expected value is {}. spearman_out.iloc[1, 2]: {}".format(
                b_v_c, spearman_out.iloc[1, 2])))

        # Should be symmetric
        self.assertEqual(spearman_out.iloc[0, 2], spearman_out.iloc[2, 0], (
            "spearman_out.iloc[0, 2]: {}".format(
                spearman_out.iloc[0, 2], spearman_out.iloc[2, 0])))

        # Diagonal should be 1s
        self.assertTrue(np.allclose(np.diagonal(spearman_out), np.ones(
            spearman_out.shape[0])), ("np.diagonal(spearman_out): {}".format(
                np.diagonal(spearman_out))))

        # Pearson
        data_df2 = pd.DataFrame(
            [[3, 7, 1, 3, 5], [9, 3, 0, 1, np.nan], [np.nan, 4, 5, 5, 1]],
            index=["a", "b", "c"], columns=["d", "e", "f", "g", "h"])
        pearson_out = steep.compute_similarity_matrix(data_df2, method="pearson")

        # This entry should be NaN
        self.assertTrue(np.isnan(pearson_out.iloc[0, 4]), (
            "Expected value is NaN. pearson_out.iloc[0, 4]: {}".format(
                pearson_out.iloc[0, 4])))
        # Should be symmetric
        self.assertEqual(pearson_out.iloc[1, 3], pearson_out.iloc[3, 1], (
            "pearson_out.iloc[1, 3]: {}".format(
                pearson_out.iloc[1, 3], pearson_out.iloc[3, 1])))
        # Diagonal should be 1s
        self.assertTrue(np.allclose(np.diagonal(pearson_out), np.ones(
            pearson_out.shape[0])), ("np.diagonal(pearson_out): {}".format(
                np.diagonal(pearson_out))))

    def test_symmetrize_if_needed(self):
        # Assymetric
        sim_df = pd.DataFrame([[0.3, 0.4, 0.1], [-0.1, 0.2, np.nan], [1, 9, 4]],
                              index=["a", "b", "c"],
                              columns=["a", "b", "c"])
        e_out_df = pd.DataFrame(
            [[0.30, 0.15, 0.55], [0.15, 0.2, np.nan], [0.55, np.nan, 4.0]])

        out_df = steep.symmetrize_if_needed(sim_df)
        self.assertTrue(np.allclose(out_df, e_out_df, equal_nan=True), (
            "out_df: \n{}".format(out_df)))

        # Symmetric
        sym_df = pd.DataFrame(
            [[0.1, 0.5, np.nan], [0.5, -0.1, 1.2], [np.nan, 1.2, 0.3]])
        out_df2 = steep.symmetrize_if_needed(sym_df)
        self.assertTrue(np.allclose(out_df2, sym_df, equal_nan=True), (
            "out_df2: \n{}".format(out_df2)))

    def test_create_ids_from_metadata(self):
        meta_df = pd.DataFrame.from_dict(
            {"cell":["a", "b", "c"], "dose":["10", "10", "1"],
             "time":["24", "24", "24"]})
        fields = ["cell", "time"]
        e_ids = ["a_24", "b_24", "c_24"]
        e_out_df = pd.DataFrame.from_dict(
            {"cell":["a", "b", "c"], "time":["24", "24", "24"]})
        e_out_df.index = pd.Index(e_ids)
        logger.debug("\ne_out_df:\n{}".format(e_out_df))

        (out_ids, out_df) = steep.create_ids_from_metadata(meta_df, "column", fields)
        self.assertEqual(out_ids, e_ids, "out_ids are incorrect: {}".format(out_ids))
        self.assertTrue(np.array_equal(out_df, e_out_df), "\nout_df:\n{}\ne_out_df:\n{}".format(
            out_df, e_out_df))

        meta_df2 = pd.DataFrame([
            ["p1", "p2", "p3"], ["-666", "-666", "-666"], ["E", "F", "G"]],
            index = ["pr_thingie", "one_more", "pr_other"])
        fields2 = np.array(["pr_thingie", "pr_other"])
        e_ids2 = ["p1_E", "p2_F", "p3_G"]
        e_out_df2 = pd.DataFrame([["p1", "p2", "p3"], ["E", "F", "G"]],
                                 index = ["pr_thingie", "pr_other"],
                                 columns = e_ids2)
        logger.debug("\ne_out_df2:\n{}".format(e_out_df2))

        (out_ids2, out_df2) = steep.create_ids_from_metadata(meta_df2, "row", fields2)
        self.assertEqual(out_ids2, e_ids2, "out_ids2 are incorrect: {}".format(out_ids2))
        self.assertTrue(np.array_equal(out_df2, e_out_df2), "\nout_df2:\n{}\ne_out_df2:\n{}".format(
            out_df2, e_out_df2))

    def test_compute_connectivity(self):
        sim_df = pd.DataFrame(
            [[1, 0.5, 1.0, -0.4, 1.1, -0.6, 0.4, -0.2, 0.1],
             [0.5, 1, 1.2, -0.8, -0.9, 0.4, 0.1, 0.1, -0.1],
             [1.0, 1.2, 1, 0.1, 0.3, 1.3, -0.9, 1.1, -0.1],
             [-0.4, -0.8, 0.1, 1, 0.5, -0.2, 0.7, 0.8, -0.3],
             [1.1, -0.9, 0.3, 0.5, 1, 0.7, 1.6, 0.4, 0.6],
             [-0.6, 0.4, 1.3, -0.2, 0.7, 1, 0.5, -0.5, -0.3],
             [0.4, 0.1, -0.9, 0.7, 1.6, 0.5, 1, 0.9, 0.6],
             [-0.2, 0.1, 1.1, 0.8, 0.4, -0.5, 0.9, 1, -0.3],
             [0.1, -0.1, -0.1, -0.3, 0.6, -0.3, 0.6, -0.3, 1]],
            index=["A", "A", "A", "B", "B", "B", "C", "C", "C"],
            columns=["A", "A", "A", "B", "B", "B", "C", "C", "C"])
        meta_df = pd.DataFrame(
            [["A", "A375", "24", "X1"], ["A", "A375", "24", "X2"],
             ["A", "A375", "24", "X3"], ["B", "PC3", "24", "X2"],
             ["B", "PC3", "24", "X2"], ["B", "PC3", "24", "X3"],
             ["C", "HT29", "24", "X1"], ["C", "HT29", "24", "X2"],
             ["C", "HT29", "24", "X3"]],
            index=["A", "A", "A", "B", "B", "B", "C", "C", "C"],
            columns=["pert", "cell", "time", "rep_num"])
        sim_gct = GCToo.GCToo(
            data_df=sim_df, row_metadata_df=meta_df, col_metadata_df=meta_df)

        e_out_df_unpivoted = pd.DataFrame.from_dict({
            "query": ["A_A375_24", "A_A375_24", "A_A375_24", "B_PC3_24",
                     "B_PC3_24", "B_PC3_24", "C_HT29_24", "C_HT29_24", "C_HT29_24"],
            "target": ["A_A375_24", "B_PC3_24", "C_HT29_24", "A_A375_24",
                      "B_PC3_24", "C_HT29_24", "A_A375_24", "B_PC3_24", "C_HT29_24"],
            "ks_statistic": [0.83, 0.36, 0.47, 0.36, 0.39, 0.33, 0.36, 0.25, 0.38],
            "p_value": [0.02, 0.42, 0.14, 0.42, 0.72, 0.52, 0.42, 0.85, 0.72],
            "ks_statistic_directed": [0.83, -0.36, -0.47, 0.36, 0.39, 0.33, -0.36, 0.25, 0.38]})
        e_out_df_unpivoted = e_out_df_unpivoted[
            ["query", "target", "ks_statistic", "p_value", "ks_statistic_directed"]]
        e_out_meta_df = pd.DataFrame(
            [["A", "A375", "24"], ["B", "PC3", "24"], ["C", "HT29", "24"]],
            index=["A_A375_24", "B_PC3_24", "C_HT29_24"],
            columns=["pert", "cell", "time"])

        (out_df_unpivoted, out_meta_df) = steep.compute_connectivity(
            sim_gct, ["pert", "cell", "time"])
        logger.debug("\nout_df_unpivoted:\n{}".format(out_df_unpivoted))

        # Check that out_meta_df is correct
        self.assertTrue(np.array_equal(e_out_meta_df, out_meta_df), (
            "\ne_out_meta_df:\n{}\nout_meta_df:\n{}").format(
            e_out_meta_df, out_meta_df))

        # Check that column names of out_df_unpivoted are correct
        self.assertEqual(list(out_df_unpivoted.columns), list(e_out_df_unpivoted.columns), (
            "list(out_df_unpivoted.columns): {}\nlist(e_out_df_unpivoted.columns): {}\n".format(
                list(out_df_unpivoted.columns), list(e_out_df_unpivoted.columns))))

        # Check that column values of out_df_unpivoted are correct
        count = 0
        self.assertTrue(np.array_equal(
            e_out_df_unpivoted.iloc[:, count], out_df_unpivoted.iloc[:, count]),
            "\ne_out_df_unpivoted[{}]:\n{}\nout_df_unpivoted[{}]:\n{}".format(
            count, e_out_df_unpivoted.iloc[:, count], count, out_df_unpivoted.iloc[:,count]))
        count = 1
        self.assertTrue(np.array_equal(
            e_out_df_unpivoted.iloc[:, count], out_df_unpivoted.iloc[:, count]),
            "\ne_out_df_unpivoted[{}]:\n{}\nout_df_unpivoted[{}]:\n{}".format(
            count, e_out_df_unpivoted.iloc[:, count], count, out_df_unpivoted.iloc[:, count]))
        count = 2
        self.assertTrue(np.allclose(
            e_out_df_unpivoted.iloc[:, count], out_df_unpivoted.iloc[:, count], atol=1e-2),
            "\ne_out_df_unpivoted[{}]:\n{}\nout_df_unpivoted[{}]:\n{}".format(
            count, e_out_df_unpivoted.iloc[:, count], count, out_df_unpivoted.iloc[:, count]))
        count = 3
        self.assertTrue(np.allclose(
            e_out_df_unpivoted.iloc[:, count], out_df_unpivoted.iloc[:, count], atol=1e-2),
            "\ne_out_df_unpivoted[{}]:\n{}\nout_df_unpivoted[{}]:\n{}".format(
            count, e_out_df_unpivoted.iloc[:, count], count, out_df_unpivoted.iloc[:, count]))
        count = 4
        self.assertTrue(np.allclose(
            e_out_df_unpivoted.iloc[:, count], out_df_unpivoted.iloc[:, count], atol=1e-1),
            "\ne_out_df_unpivoted[{}]:\n{}\nout_df_unpivoted[{}]:\n{}".format(
            count, e_out_df_unpivoted.iloc[:, count], count, out_df_unpivoted.iloc[:, count]))

    def test_create_distributions_for_ks_test(self):
        # 2 perts; 3 replicates each
        in_df = pd.DataFrame(
            [[np.nan, 0.5, 1.0, -0.4, 1.1, -0.6],
             [0.5, np.nan, 1.2, -0.8, -0.9, 0.4],
             [1.0, 1.2, np.nan, 0.1, 0.3, 1.3],
             [-0.4, -0.8, 0.1, np.nan, 0.5, -0.2],
             [1.1, -0.9, 0.3, 0.5, np.nan, 0.7],
             [-0.6, 0.4, 1.3, -0.2, 0.7, np.nan]],
            index=["A", "A", "A", "B", "B", "B"],
            columns=["A", "A", "A", "B", "B", "B"])
        logger.debug("\nin_df:\n{}".format(in_df))
        e_queries = ["A", "A", "B", "B"]
        e_targets = ["A", "B", "A", "B"]
        e_tests = [[0.5, 1.0, 1.2],
                   [-0.4, -0.8, 0.1, 1.1, -0.9, 0.3, -0.6, 0.4, 1.3],
                   [-0.4, -0.8, 0.1, 1.1, -0.9, 0.3, -0.6, 0.4, 1.3],
                   [0.5, -0.2, 0.7]]
        e_nulls = [[-0.4, -0.8, 0.1, 1.1, -0.9, 0.3, -0.6, 0.4, 1.3],
                   [0.5, -0.2, 0.7], [0.5, 1.0, 1.2],
                   [-0.4, -0.8, 0.1, 1.1, -0.9, 0.3, -0.6, 0.4, 1.3]]
        e_unique_perts = ["A", "B"]

        [queries, targets, tests, nulls, unique_perts] = steep.create_distributions_for_ks_test(in_df)

        self.assertItemsEqual(e_queries, queries, "\ne_queries:\n{}\nqueries:\n{}".format(
            e_queries, queries))
        self.assertItemsEqual(e_targets, targets, "\ne_targets:\n{}\ntargets:\n{}".format(
            e_targets, targets))
        self.assertItemsEqual(e_unique_perts, unique_perts, "\ne_unique_perts:\n{}\nunique_perts:\n{}".format(
            e_unique_perts, unique_perts))
        for count, _ in enumerate(e_tests):
            self.assertItemsEqual(e_tests[count], tests[count], (
                "\ne_tests[{}]:\n{}\ntests[{}]:\n{}").format(
                count, e_tests[count], count, tests[count]))
        for count, _ in enumerate(e_nulls):
            self.assertItemsEqual(e_nulls[count], nulls[count], (
                "\ne_nulls[{}]:\n{}\nnulls[{}]:\n{}").format(
                count, e_nulls[count], count, nulls[count]))

        # 3 perts; 2, 2, and 1 replicates
        in_df2 = pd.DataFrame(
            [[np.nan, 0.5, 1.0, -0.4, 1.1],
             [0.5, np.nan, 1.2, -0.8, -0.9],
             [1.0, 1.2, np.nan, 0.1, 0.3],
             [-0.4, -0.8, 0.1, np.nan, 0.5],
             [1.1, -0.9, 0.3, 0.5, np.nan]],
            index=["A", "A", "B", "B", "C"],
            columns=["A", "A", "B", "B", "C"])
        logger.debug("\nin_df2:\n{}".format(in_df2))
        e_queries2 = ["A", "A", "A", "B", "B", "B", "C", "C", "C"]
        e_targets2 = ["A", "B", "C", "A", "B", "C", "A", "B", "C"]
        e_tests2 = [[], [1.0, 1.2, -0.4, -0.8], [],
                    [1.0, -0.4, 1.2, -0.8], [], [], [1.1, -0.9],
                    [0.3, 0.5], []]
        e_nulls2 = [[], [0.1, 0.3, 0.5], [], [0.5, 1.1, -0.9],
                    [], [], [0.5, 1.0, -0.4, 1.2, -0.8],
                    [1.0, 1.2, 0.1, -0.4, -0.8], []]
        e_unique_perts2 = ["A", "B", "C"]

        [queries2, targets2, tests2, nulls2, unique_perts2] = steep.create_distributions_for_ks_test(in_df2)

        self.assertItemsEqual(e_queries2, queries2, "\ne_queries2:\n{}\nqueries2:\n{}".format(
            e_queries2, queries2))
        self.assertItemsEqual(e_targets2, targets2, "\ne_targets2:\n{}\ntargets2:\n{}".format(
            e_targets2, targets2))
        self.assertItemsEqual(e_unique_perts2, unique_perts2, "\ne_unique_perts2:\n{}\nunique_perts2:\n{}".format(
            e_unique_perts2, unique_perts2))
        for count, test in enumerate(e_tests2):
            self.assertItemsEqual(e_tests2[count], tests2[count], (
                "\ne_tests2[{}]:\n{}\ntests2[{}]:\n{}").format(
                count, e_tests2[count], count, tests2[count]))
        for count, null in enumerate(e_nulls2):
            self.assertItemsEqual(e_nulls2[count], nulls2[count], (
                "\ne_nulls2[{}]:\n{}\nnulls2[{}]:\n{}").format(
                count, e_nulls2[count], count, nulls2[count]))

        # 3 perts; 3 replicates each
        in_df3 = pd.DataFrame(
            [[np.nan, 0.5, 1.0, -0.4, 1.1, -0.6, 0.4, -0.2, 0.1],
             [0.5, np.nan, 1.2, -0.8, -0.9, 0.4, 0.1, 0.1, -0.1],
             [1.0, 1.2, np.nan, 0.1, 0.3, 1.3, -0.9, 1.1, -0.1],
             [-0.4, -0.8, 0.1, np.nan, 0.5, -0.2, 0.7, 0.8, -0.3],
             [1.1, -0.9, 0.3, 0.5, np.nan, 0.7, 1.6, 0.4, 0.6],
             [-0.6, 0.4, 1.3, -0.2, 0.7, np.nan, 0.5, -0.5, -0.3],
             [0.4, 0.1, -0.9, 0.7, 1.6, 0.5, np.nan, 0.9, 0.6],
             [-0.2, 0.1, 1.1, 0.8, 0.4, -0.5, 0.9, np.nan, -0.3],
             [0.1, -0.1, -0.1, -0.3, 0.6, -0.3, 0.6, -0.3, np.nan]],
            index=["A", "A", "A", "B", "B", "B", "C", "C", "C"],
            columns=["A", "A", "A", "B", "B", "B", "C", "C", "C"])
        logger.debug("\nin_df3:\n{}".format(in_df3))
        e_queries3 = ["A", "A", "A", "B", "B", "B", "C", "C", "C"]
        e_targets3 = ["A", "B", "C", "A", "B", "C", "A", "B", "C"]
        e_tests3_0 = [0.5, 1.0, 1.2]
        e_nulls3_0 = [-0.4, 1.1, -0.6, 0.4, -0.2, 0.1, -0.8, -0.9, 0.4, 0.1,
                     0.1, -0.1, 0.1, 0.3, 1.3, -0.9, 1.1, -0.1]
        e_tests3_1 = [-0.4,-0.8, 0.1, 1.1, -0.9, 0.3, -0.6, 0.4, 1.3]
        e_nulls3_1 = [0.5, -0.2, 0.7, 0.8, -0.3, 0.7, 1.6, 0.4, 0.6, 0.5, -0.5, -0.3]
        e_unique_perts3 = ["A", "B", "C"]

        [queries3, targets3, tests3, nulls3, unique_perts3] = steep.create_distributions_for_ks_test(in_df3)

        self.assertItemsEqual(e_queries3, queries3, "\ne_queries3:\n{}\nqueries3:\n{}".format(
            e_queries3, queries3))
        self.assertItemsEqual(e_targets3, targets3, "\ne_targets3:\n{}\ntargets3:\n{}".format(
            e_targets3, targets3))
        self.assertItemsEqual(e_unique_perts3, unique_perts3, "\ne_unique_perts3:\n{}\nunique_perts3:\n{}".format(
            e_unique_perts3, unique_perts3))

        # Too big to test exhaustively
        self.assertItemsEqual(e_tests3_0, tests3[0], (
                "\ne_tests_0:\n{}\ntests3[0]:\n{}").format(
                e_tests3_0, tests3[0]))
        self.assertItemsEqual(e_tests3_1, tests3[1], (
                "\ne_tests_1:\n{}\ntests3[1]:\n{}").format(
                e_tests3_1, tests3[1]))
        self.assertItemsEqual(e_nulls3_0, nulls3[0], (
                "\ne_nulls3_0:\n{}\nnulls3[0]:\n{}").format(
                e_nulls3_0, nulls3[0]))
        self.assertItemsEqual(e_nulls3_1, nulls3[1], (
                "\ne_nulls3_1:\n{}\nnulls3[1]:\n{}").format(
                e_nulls3_1, nulls3[1]))

        # 1 pert; 3 replicates

        # Assymetric similarity matrix...


    def test_configure_out_gct_name(self):
        # Name from args is None
        in_gct = "a/b/c/my_gct.gct"
        sim_suffix = ".steep.similarity.gct"
        args_sim_name = None
        out_sim_name = steep.configure_out_gct_name(in_gct, sim_suffix, args_sim_name)
        self.assertEqual(out_sim_name, "my_gct.gct.steep.similarity.gct")

        # Name is specified in args
        args_sim_name2 = "some_other_thing.gct"
        out_sim_name2 = steep.configure_out_gct_name(
            in_gct, sim_suffix, args_sim_name2)
        self.assertEqual(out_sim_name2, "some_other_thing.gct")


    def test_ks_test(self):
        # Verify that doubling the number of elements going into the
        # distributions for a 2-sample KS-test DOES affect results
        null1 = [6, 3, 5, 2, 5]
        test1 = [4, 2, 3]
        null2 = null1 + null1
        test2 = test1 + test1

        (ks_stat, pval) = stats.ks_2samp(test1, null1)
        (ks_stat2, pval2) = stats.ks_2samp(test2, null2)

        self.assertAlmostEqual(ks_stat, ks_stat2)
        self.assertNotEqual(pval, pval2)

if __name__ == "__main__":
    setup_logger.setup(verbose=False)
    unittest.main()