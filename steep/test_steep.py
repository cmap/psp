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

# Functional tests dir lives within the dry directory
FUNCTIONAL_TESTS_DIR = "steep/functional_tests"

# Debugging purposes
import in_out.write_gctoo as wg


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

    def test_compute_similarity_against_R_code(self):
        INPUT_GCT_PATH = os.path.join(FUNCTIONAL_TESTS_DIR, "PC3_QCNORM_n189x59.gct")
        JJ_PC3_SIM_PATH = os.path.join(FUNCTIONAL_TESTS_DIR, "R_PC3_SIM_n189x189.csv")

        # Read in PC3 data and JJ's R version of the similarity matrix
        pc3_qcnorm_gct = parse_gctoo.parse(INPUT_GCT_PATH)
        jj_pc3_sim_df = pd.read_csv(JJ_PC3_SIM_PATH, index_col=0)

        # Compute similarity
        pc3_sim_df = steep.compute_similarity_matrix(pc3_qcnorm_gct.data_df, method="spearman")

        self.assertTrue(np.allclose(pc3_sim_df, jj_pc3_sim_df), (
            "\npc3_sim_df:\n{}\njj_pc3_sim_df:\n{}".format(pc3_sim_df, jj_pc3_sim_df)))

    def test_create_group_ids(self):
        # Create group_id column
        meta_df = pd.DataFrame.from_dict(
            {"cell":["a", "b", "c"], "dose":["10", "10", "1"],
             "time":["24", "24", "24"]})
        group_fields = ["cell", "time"]
        e_group_ids = ["a_24", "b_24", "c_24"]
        e_out_df = pd.DataFrame.from_dict(
            {"cell":["a", "b", "c"], "time":["24", "24", "24"],
              "group_id":e_group_ids})
        e_out_df = e_out_df[["cell", "time", "group_id"]]
        logger.debug("\ne_out_df:\n{}".format(e_out_df))

        (out_group_ids, out_df) = steep.create_group_ids(meta_df, "col", group_fields)
        self.assertEqual(out_group_ids, e_group_ids, (
            "out_group_ids: {}\ne_group_ids: {}".format(out_group_ids, e_group_ids)))
        self.assertTrue(np.array_equal(out_df, e_out_df), "\nout_df:\n{}\ne_out_df:\n{}".format(
            out_df, e_out_df))

        # Create group_id row
        meta_df2 = pd.DataFrame([
            ["p1", "p2", "p3"], ["-666", "-666", "-666"], ["E", "F", "G"]],
            index = ["pr_thingie", "one_more", "pr_other"])
        group_fields2 = np.array(["pr_thingie", "pr_other"])
        e_group_ids2 = ["p1_E", "p2_F", "p3_G"]
        e_out_df2 = pd.DataFrame([
            ["p1", "p2", "p3"], ["E", "F", "G"], e_group_ids2],
            index = ["pr_thingie", "pr_other", "group_id"])
        logger.debug("\ne_out_df2:\n{}".format(e_out_df2))

        (out_group_ids2, out_df2) = steep.create_group_ids(meta_df2, "row", group_fields2)
        self.assertEqual(out_group_ids2, e_group_ids2, (
            "out_group_ids2: {}\ne_group_ids2: {}".format(out_group_ids2, e_group_ids2)))
        self.assertTrue(np.array_equal(out_df2, e_out_df2), "\nout_df2:\n{}\ne_out_df2:\n{}".format(
            out_df2, e_out_df2))

    def test_compute_connectivity(self):
        sim_df = pd.DataFrame(
            [[1, 0.5, 1.0, -0.4, 1.1, -0.6, 0.4, -0.2, 0.1, 0.5, -0.2, 0.1],
             [0.5, 1, 1.2, -0.8, -0.9, 0.4, 0.1, 0.1, -0.1, -0.7, 0.3, 0.4],
             [1.0, 1.2, 1, 0.1, 0.3, 1.3, -0.9, 1.1, -0.1, -0.2, 0.2, 0.2],
             [-0.4, -0.8, 0.1, 1, 0.5, -0.2, 0.7, 0.8, -0.3, -1, 0.1, 0.6],
             [1.1, -0.9, 0.3, 0.5, 1, 0.7, 1.6, 0.4, 0.6, 0.4, 0.4, 0.4],
             [-0.6, 0.4, 1.3, -0.2, 0.7, 1, 0.5, -0.5, -0.3, 0.2, -0.9, -0.1],
             [0.4, 0.1, -0.9, 0.7, 1.6, 0.5, 1, 0.9, 0.6, 0.6, -0.5, -0.1],
             [-0.2, 0.1, 1.1, 0.8, 0.4, -0.5, 0.9, 1, -0.3, 0.6, -0.7, 0.8],
             [0.1, -0.1, -0.1, -0.3, 0.6, -0.3, 0.6, -0.3, 1, 0.2, 0.8, -0.9],
             [0.5, -0.7, -0.2, -1, 0.4, 0.2, 0.6, 0.6, 0.2, 1, -0.2, 0.3],
             [-0.2, 0.3, 0.2, 0.1, 0.4, -0.9, -0.5, -0.7, 0.8, -0.2, 1, -0.8],
             [0.1, 0.4, 0.2, 0.6, 0.4, -0.1, -0.1, 0.8, -0.9, 0.3, -0.8, 1]],
            index=["a1", "a2", "a3", "b1", "b2", "b3",
                   "c1", "c2", "c3", "d1", "d2", "d3"],
            columns=["a1", "a2", "a3", "b1", "b2", "b3",
                   "c1", "c2", "c3", "d1", "d2", "d3"])
        meta_df = pd.DataFrame(
            [["DMSO", "PC3", "24", "X1"], ["DMSO", "PC3", "24", "X2"],
             ["DMSO", "PC3", "24", "X3"], ["DMSO", "MCF7", "24", "X1"],
             ["DMSO", "MCF7", "24", "X2"], ["DMSO", "MCF7", "24", "X3"],
             ["V", "PC3", "24", "X1"], ["V", "PC3", "24", "X2"],
             ["V", "PC3", "24", "X3"], ["V", "MCF7", "24", "X1"],
             ["V", "MCF7", "24", "X2"], ["V", "MCF7", "24", "X3"]],
            index=["a1", "a2", "a3", "b1", "b2", "b3",
                   "c1", "c2", "c3", "d1", "d2", "d3"],
            columns=["pert_iname", "cell_id", "pert_time", "rep_num"])
        sim_gct = GCToo.GCToo(
            data_df=sim_df, row_metadata_df=meta_df, col_metadata_df=meta_df)

        # Debugging purposes
        wg.write(sim_gct, "small_sim_gct_for_testing_connectivity.gct")

        e_out_df_unpivoted = pd.DataFrame.from_dict({
            "query": ["DMSO_MCF7_24", "DMSO_MCF7_24", "DMSO_MCF7_24", "DMSO_MCF7_24",
                      "DMSO_PC3_24", "DMSO_PC3_24", "DMSO_PC3_24", "DMSO_PC3_24",
                      "V_MCF7_24", "V_MCF7_24", "V_MCF7_24", "V_MCF7_24",
                      "V_PC3_24", "V_PC3_24", "V_PC3_24", "V_PC3_24"],
            "target": ["DMSO_MCF7_24", "DMSO_PC3_24", "V_MCF7_24", "V_PC3_24",
                       "DMSO_MCF7_24", "DMSO_PC3_24", "V_MCF7_24", "V_PC3_24",
                       "DMSO_MCF7_24", "DMSO_PC3_24", "V_MCF7_24", "V_PC3_24",
                       "DMSO_MCF7_24", "DMSO_PC3_24", "V_MCF7_24", "V_PC3_24"],
            "ks_statistic": [0.37, 0.35, 0.17, 0.29, 0.30, 0.85, 0.23, 0.41, 0.31, 0.24, 0.37, 0.24, 0.32, 0.35, 0.40, 0.33],
            "p_value": [0.75, 0.35, 0.98, 0.60, 0.53, 0.02,0.81, 0.18, 0.47, 0.81, 0.75, 0.81, 0.47, 0.35, 0.21, 0.85],
            "ks_statistic_signed": [0.37, 0.35, 0.17, 0.29, -0.30, 0.85, 0.23, -0.41, -0.31, 0.24, -0.37, 0.24, 0.32, -0.35, 0.40, 0.33]})
        e_out_df_unpivoted = e_out_df_unpivoted[
            ["query", "target", "ks_statistic", "p_value", "ks_statistic_signed"]]
        e_out_meta_df = pd.DataFrame(
            [["DMSO", "PC3", "24"], ["DMSO", "MCF7", "24"],
             ["V", "PC3", "24"], ["V", "MCF7", "24"]],
            index=["DMSO_MCF7_24", "DMSO_PC3_24", "V_MCF7_24", "V_PC3_24"],
            columns=["pert_iname", "cell_id", "pert_time"])

        (out_df_unpivoted, out_meta_df) = steep.compute_connectivity(
            sim_gct, ["pert_iname", "cell_id", "pert_time"])
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
        in_df4 = pd.DataFrame(
            [[np.nan, 0.5, 1.0],
             [0.5, np.nan, 1.2],
             [1.0, 1.2, np.nan]],
            index=["A", "A", "A"],
            columns=["A", "A", "A"])
        logger.debug("\nin_df4:\n{}".format(in_df4))
        e_queries4 = ["A"]
        e_targets4 = ["A"]
        e_tests4 = [[]]
        e_nulls4 = [[]]
        e_unique_perts4 = ["A"]
        
        [queries4, targets4, tests4, nulls4, unique_perts4] = steep.create_distributions_for_ks_test(in_df4)
        
        self.assertItemsEqual(e_queries4, queries4, "\ne_queries4:\n{}\nqueries4:\n{}".format(
            e_queries4, queries4))
        self.assertItemsEqual(e_targets4, targets4, "\ne_targets4:\n{}\ntargets4:\n{}".format(
            e_targets4, targets4))
        self.assertItemsEqual(e_tests4, tests4, "\ne_tests4:\n{}\ntests4:\n{}".format(
            e_tests4, tests4))
        self.assertItemsEqual(e_nulls4, nulls4, "\ne_nulls4:\n{}\nnulls4:\n{}".format(
            e_nulls4, nulls4))
        self.assertItemsEqual(e_unique_perts4, unique_perts4, "\ne_unique_perts4:\n{}\nunique_perts4:\n{}".format(
            e_unique_perts4, unique_perts4))

    def test_assemble_output_conn_gcts(self):
        conn_df_unpivoted = pd.DataFrame([
            ["k", "J", 0.1, 0.03, 0.1],
            ["J", "k", 0.2, 0.04, -0.2],
            ["k", "l", 0.3, 0.05, 0.3],
            ["k", "k", 0.4, 0.04, -0.4],
            ["l", "l", 0.5, 0.03, 0.5],
            ["J", "J", 0.6, 0.03, 0.6],
            ["l", "k", 0.8, 0.2, 0.8],
            ["l", "J", 0.7, 0.3, 0.7],
            ["J", "l", 0.9, 0.1, -0.9]],
            columns=["query", "target", "ks_statistic", "p_value", "ks_statistic_signed"])
        conn_meta_df = pd.DataFrame([["l1", "l2"], ["j1", "j2"], ["k1", "k2"]],
            index=["l", "J", "k"], columns=["field1", "field2"])
        e_conn_data_df = pd.DataFrame(
            [[0.6, 0.2, 0.9], [0.1, 0.4, 0.3], [0.7, 0.8, 0.5]],
            index=["J", "k", "l"], columns=["J", "k", "l"])
        e_pval_data_df = pd.DataFrame(
            [[0.03, 0.04, 0.1], [0.03, 0.04, 0.05], [0.3, 0.2, 0.03]],
            index=["J", "k", "l"], columns=["J", "k", "l"])
        e_conn_signed_data_df = pd.DataFrame(
            [[0.6, -0.2, -0.9], [0.1, -0.4, 0.3], [0.7, 0.8, 0.5]],
            index=["J", "k", "l"], columns=["J", "k", "l"])
        e_conn_meta_df = pd.DataFrame([["j1", "j2"], ["k1", "k2"], ["l1", "l2"]],
            index=["J", "k", "l"], columns=["field1", "field2"])

        (conn_gct, pval_gct, signed_conn_gct) = steep.assemble_output_conn_gcts(
            conn_df_unpivoted, conn_meta_df)

        # Check matrices
        self.assertTrue(conn_gct.data_df.equals(e_conn_data_df), (
            "\nconn_gct.data_df:\n{}\ne_conn_data_df:\n{}".format(
                conn_gct.data_df, e_conn_data_df)))
        self.assertTrue(pval_gct.data_df.equals(e_pval_data_df), (
            "\npval_gct.data_df:\n{}\ne_pval_data_df:\n{}".format(
                pval_gct.data_df, e_pval_data_df)))
        self.assertTrue(signed_conn_gct.data_df.equals(e_conn_signed_data_df), (
            "\nsigned_conn_gct.data_df:\n{}\ne_conn_signed_data_df:\n{}".format(
                signed_conn_gct.data_df, e_conn_signed_data_df)))

        # Check metadata
        self.assertTrue(conn_gct.row_metadata_df.equals(e_conn_meta_df))
        self.assertTrue(conn_gct.col_metadata_df.equals(e_conn_meta_df))
        self.assertTrue(pval_gct.row_metadata_df.equals(e_conn_meta_df))
        self.assertTrue(pval_gct.col_metadata_df.equals(e_conn_meta_df))
        self.assertTrue(signed_conn_gct.row_metadata_df.equals(e_conn_meta_df))
        self.assertTrue(signed_conn_gct.col_metadata_df.equals(e_conn_meta_df))

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