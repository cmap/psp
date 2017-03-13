import logging
import unittest
import igraph as ig
import pandas as pd
import tasseography

import broadinstitute_cmap.io.pandasGEXpress.GCToo as GCToo
import broadinstitute_psp.utils.setup_logger as setup_logger

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestTasseography(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.sym_data_df = pd.DataFrame(
            [[0.9, 0.4, 0.6], [0.4, 1.0, -0.3], [0.6, -0.3, 1.0]],
            index=["a", "b", "c"], columns=["a", "b", "c"])
        cls.sym_meta_df = pd.DataFrame(
            [["A375", "great"], ["A375", "bad"], ["A375", "ok"]],
            index=["a", "b", "c"], columns=["cell_id", "pert_type"])
        cls.sym_gct = GCToo.GCToo(cls.sym_data_df, cls.sym_meta_df, cls.sym_meta_df)

        cls.asym_data_df = pd.DataFrame(
            [[0.1, 0.4], [0.7, -0.1], [0.9, 0.9]],
            index=["A", "B", "C"], columns=["o", "k"])
        cls.asym_row_meta_df = pd.DataFrame(
            ["3h", "1h", "2h"], index=["A", "B", "C"], columns=["pert_time"])
        cls.asym_col_meta_df = pd.DataFrame(
            ["MCF7", "A375"], index=["o", "k"], columns=["cell_id"])
        cls.asym_gct = GCToo.GCToo(cls.asym_data_df, cls.asym_row_meta_df, cls.asym_col_meta_df)

        cls.sym_g = ig.Graph()
        cls.sym_g.add_vertices(3)
        cls.sym_g.add_edges([(0, 1), (0, 2), (1, 2)])
        cls.sym_g.vs["id"] = ["a", "b", "c"]
        cls.sym_g.vs["cell_id"] = ["A375", "A375", "A375"]
        cls.sym_g.vs["pert_type"] = ["great", "bad", "ok"]
        cls.sym_g.es["weight"] = [0.4, 0.6, -0.3]

        cls.asym_g = ig.Graph()
        cls.asym_g.add_vertices(5)
        cls.asym_g.add_edges([(0, 3), (0, 4), (1, 3), (1, 4), (2, 3), (2, 4)])
        cls.asym_g.vs["id"] = ["A", "B", "C", "o", "k"]
        cls.asym_g.vs["cell_id"] = [None, None, None, "MCF7", "A375"]
        cls.asym_g.vs["pert_time"] = ["3h", "1h", "2h", None, None]
        cls.asym_g.es["weight"] = [0.1, 0.4, 0.7, -0.1, 0.9, 0.9]

    def test_sym_gct_to_graph(self):
        logger.debug("self.sym_gct:\n{}".format(self.sym_gct))

        with self.assertRaises(Exception) as e:
            tasseography.sym_gct_to_graph(self.asym_gct, [])
        self.assertIn("Row metadata must", str(e.exception))

        out = tasseography.sym_gct_to_graph(self.sym_gct, ["pert_type", "cell_id"])
        logger.debug("out: {}".format(out))

        self.assertItemsEqual(out.vs["id"], self.sym_g.vs["id"])
        self.assertItemsEqual(out.vs["cell_id"], self.sym_g.vs["cell_id"])
        self.assertItemsEqual(out.vs["pert_type"], self.sym_g.vs["pert_type"])
        self.assertItemsEqual(out.es["weight"], self.sym_g.es["weight"])

    def test_asym_gct_to_graph(self):
        logger.debug("self.asym_gct:\n{}".format(self.asym_gct))

        with self.assertRaises(Exception) as e:
            tasseography.asym_gct_to_graph(self.asym_gct, ["pert_time"], ["pert_time"])
        self.assertIn("field pert_time not in column metadata", str(e.exception))

        out = tasseography.asym_gct_to_graph(self.asym_gct, ["pert_time"], ["cell_id"])
        logger.debug("out: {}".format(out))

        self.assertItemsEqual(out.vs["id"], self.asym_g.vs["id"])
        self.assertItemsEqual(out.vs["cell_id"], self.asym_g.vs["cell_id"])
        self.assertItemsEqual(out.vs["pert_time"], self.asym_g.vs["pert_time"])
        self.assertItemsEqual(out.es["weight"], self.asym_g.es["weight"])

    def test_melt_df(self):
        expected = pd.DataFrame(
            [["A", "o", 0.1], ["B", "o", 0.7], ["C", "o", 0.9],
             ["A", "k", 0.4], ["B", "k", -0.1], ["C", "k", 0.9]],
            columns=["row", "column", "value"])

        out = tasseography.melt_df(self.asym_data_df)
        pd.util.testing.assert_frame_equal(expected, out)


    def test_remove_edges_below_thresh(g, thresh):
        pass

    def test_extract_subgraph(g, node_field, node_value):
        pass

    def test_plot_network(g, layout, node_color_field, edge_color_field):
        pass

    def test_plot_bipartite(g):
        pass

    def test_graph_to_text_files(g, out_node_name, out_edge_name):
        pass


if __name__ == '__main__':
    setup_logger.setup(verbose=True)
    unittest.main()
