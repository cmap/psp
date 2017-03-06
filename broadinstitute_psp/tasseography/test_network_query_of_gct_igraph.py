import igraph as ig
import logging
import pandas as pd
import os
import unittest
import network_query_of_gct_igraph as nq
import broadinstitute_psp.utils.setup_logger as setup_logger

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestNetworkQueryOfGctIgraph(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.g = ig.Graph()
        cls.g.add_vertices(4)
        cls.g.add_edges([(0, 2), (1, 2)])
        cls.g.vs["id"] = ["a", "b", "c", "d"]
        cls.g.vs["moa"] = ["J", "K", "L", "M"]
        cls.g.es["weight"] = [0.8, 0.85]

    def test_make_color_dict(self):
        meta_df = pd.DataFrame(
            [["A375", "A"],
             ["A375", "B"],
             ["A375", "A"]],
            index=["a", "b", "c"],
            columns=["cell", "moa"])
        e_color_dict = {'A': (1.0, 0.0, 0.0, 1.0),  # red
                        'B': (0.0, 1.0, 0.0, 1.0)}  # blue
        color_dict = nq.make_color_dict(meta_df, "moa")
        self.assertEqual(color_dict, e_color_dict)

    def test_write_edge_tsv(self):
        e_edge_df = pd.DataFrame([["a", "c", 0.8], ["b", "c", 0.85]],
                                 columns=["query", "target", "value"])
        edge_name = "test_network_query_out_edge.tsv"
        nq.write_edge_tsv(self.g, edge_name)
        edge_df = pd.read_csv(edge_name, sep="\t")
        pd.util.testing.assert_frame_equal(edge_df, e_edge_df)

        os.remove(edge_name)

    def test_write_node_tsv(self):
        e_node_df = pd.DataFrame([["a", "J"], ["b", "K"],
                                  ["c", "L"], ["d", "M"]],
                                 columns=["id", "moa"])
        node_name = "test_network_query_out_node.tsv"
        nq.write_node_tsv(self.g, node_name)
        node_df = pd.read_csv(node_name, sep="\t")
        pd.util.testing.assert_frame_equal(node_df, e_node_df)

        os.remove(node_name)


if __name__ == '__main__':
    setup_logger.setup(verbose=True)
    unittest.main()
