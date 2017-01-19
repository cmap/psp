import logging
import pandas as pd
import unittest
import broadinstitute_psp.utils.setup_logger as setup_logger
import network_query_of_gct as nq

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestNetworkQuery(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.melted_df = pd.DataFrame(
            [["a", "b", 0.8],
             ["a", "c", 0.4],
             ["b", "a", 0.65],
             ["b", "c", 0.55],
             ["c", "a", 0.4],
             ["c", "b", 0.9]],
            index=range(6), columns=["query", "target", "value"])
        cls.col_meta = pd.DataFrame(
            [["a_moa", "a_name"],
             ["b_moa", "b_name"],
             ["c_moa", "c_name"]],
            index=["a", "b", "c"], columns=["MoA", "pert_iname"])

    def test_get_all_cxns_above_thresh(self):
        e_edge_df = pd.DataFrame(
            [["a", "b", 0.8, 0.8, 1.], ["c", "b", 0.9, 0.9, 1.]],
            index=[0, 5], columns=["query", "target", "value", "abs_value", "sign"])
        e_node_df = self.col_meta.copy()
        e_node_df.index.name = "pert_id"
        e_node_df = e_node_df.loc[["a", "c", "b"], ["pert_iname", "MoA"]]

        edge_df, node_df = nq.get_all_cxns_above_thresh(
            self.melted_df, self.col_meta, 0.8, ["pert_iname", "MoA"])

        pd.util.testing.assert_frame_equal(edge_df, e_edge_df)
        pd.util.testing.assert_frame_equal(node_df, e_node_df)

    def test_convert_pert_iname_to_pert_id(self):
        e_pert_id = "b"
        pert_id = nq.convert_pert_iname_to_pert_id("b_name", self.col_meta)
        self.assertEqual(pert_id, e_pert_id)

        with self.assertRaises(AssertionError) as e:
            nq.convert_pert_iname_to_pert_id("b_nae", self.col_meta)
        self.assertIn("b_nae not found in column", str(e.exception))

    def test_create_edge_df(self):
        edge_df1 = pd.DataFrame(
            [["a", "b", 0.8], ["c", "b", 0.9]],
            index=[0, 5], columns=["query", "target", "value"])
        edge_df2 = pd.DataFrame(
            [["d", "e", -0.8], ["a", "b", 0.8]],
            index=[6, 0], columns=["query", "target", "value"])
        e_edge_df = pd.DataFrame(
            [["a", "b", 0.8, 0.8, 1.0], ["c", "b", 0.9, 0.9, 1.0], ["d", "e", -0.8, 0.8, -1.0]],
            index=[0, 5, 6], columns=["query", "target", "value", "abs_value", "sign"])

        edge_df = nq.create_edge_df([edge_df1, edge_df2])
        pd.util.testing.assert_frame_equal(edge_df, e_edge_df)

    def test_get_cxns_to_query_above_thresh(self):
        e_cxns = self.melted_df.iloc[[0], :]
        e_targets = ["b"]
        cxns, targets = nq.get_cxns_to_query_above_thresh(self.melted_df, "a", 0.7)

        pd.util.testing.assert_frame_equal(cxns, e_cxns)
        self.assertEqual(targets, e_targets)

    def test_get_second_order_cxns(self):
        first_order_targets = ["b", "c"]
        e_cxns = self.melted_df.iloc[[5], :]
        cxns = nq.get_second_order_cxns(self.melted_df, first_order_targets, 0.7)

        pd.util.testing.assert_frame_equal(cxns, e_cxns)

    def test_network_query_of_df(self):
        more_complex_melted_df = pd.DataFrame(
            [["a", "b", 0.8],
             ["a", "c", 0.4],
             ["a", "d", -0.65],
             ["b", "a", 0.5],
             ["b", "c", 0.55],
             ["b", "d", 0.95],
             ["c", "a", 0.4],
             ["c", "b", 0.9],
             ["c", "d", 0.9],
             ["d", "a", -0.7],
             ["d", "b", 0.6],
             ["d", "c", 0.1]],
            index=range(12), columns=["query", "target", "value"])
        more_complex_col_meta = pd.DataFrame(
            [["a_moa", "a_name"],
             ["d_moa", "d_name"],
             ["b_moa", "b_name"],
             ["c_moa", "c_name"]],
            index=["a", "d", "b", "c"], columns=["moa", "pert_iname"])

        e_edge_df = more_complex_melted_df.iloc[[9, 10, 0], :]
        e_edge_df["abs_value"] = [0.7, 0.6, 0.8]
        e_edge_df["sign"] = [-1., 1., 1.]
        e_node_df = more_complex_col_meta.loc[["a", "b", "d"]]
        e_node_df = e_node_df[["pert_iname", "moa"]]
        e_node_df.index.name = "pert_id"

        node_df, edge_df = nq.network_query_of_df(
            more_complex_melted_df, more_complex_col_meta, ["d_name"], 0.6, ["pert_iname", "moa"])

        pd.util.testing.assert_frame_equal(edge_df, e_edge_df)
        pd.util.testing.assert_frame_equal(node_df, e_node_df)


if __name__ == '__main__':
    setup_logger.setup(verbose=True)
    unittest.main()
