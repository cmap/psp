import logging
import unittest
import igraph as ig
import numpy as np
import os
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
            [[0.1, 0.4], [-0.7, -0.1], [np.nan, 0.9]],
            index=["A", "B", "C"], columns=["o", "k"])
        cls.asym_row_meta_df = pd.DataFrame(
            ["3h", "1h", "2h"], index=["A", "B", "C"], columns=["pert_time"])
        cls.asym_col_meta_df = pd.DataFrame(
            [["MCF7", "3h"], ["A375", "6h"]], index=["o", "k"],
            columns=["cell_id", "pert_time"])
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
        cls.asym_g.vs["type"] = [False, False, False, True, True]
        cls.asym_g.vs["cell_id"] = [None, None, None, "MCF7", "A375"]
        cls.asym_g.vs["pert_time"] = ["3h", "1h", "2h", "3h", "6h"]
        cls.asym_g.es["weight"] = [0.1, 0.4, -0.7, -0.1, np.nan, 0.9]

    def test_main_sym(self):
        out_fig_name = "test_main_sym_output.png"
        out_gml_name = "test_main_sym_output.gml"

        tasseography.main_sym(
            self.sym_gct, out_fig_name, out_gml_name, ["cell_id", "pert_type"],
            ["great"], "pert_type", 0, "pert_type", None, layout="fr")

        # Make sure plot and gml files were produced
        self.assertTrue(os.path.exists(out_fig_name))
        self.assertTrue(os.path.exists(out_gml_name))

        # Now, remove them
        os.remove(out_fig_name)
        os.remove(out_gml_name)

    def test_main_asym(self):
        out_fig_name = "test_main_asym_output.png"
        out_gml_name = "test_main_asym_output.gml"

        tasseography.main_asym(
            self.asym_gct, out_fig_name, out_gml_name, ["pert_time"],
            ["cell_id", "pert_time"], ["1h"], "pert_time", "row", 0,
            "pert_time", None)

        # Make sure plot and gml files were produced
        self.assertTrue(os.path.exists(out_fig_name))
        self.assertTrue(os.path.exists(out_gml_name))

        # Now, remove them
        os.remove(out_fig_name)
        os.remove(out_gml_name)

    def test_sym_gct_to_graph(self):
        logger.debug("self.sym_gct:\n{}".format(self.sym_gct))

        with self.assertRaises(Exception) as e:
            tasseography.sym_gct_to_graph(self.asym_gct, [])
        self.assertIn("Row metadata must", str(e.exception))

        out = tasseography.sym_gct_to_graph(self.sym_gct, ["pert_type", "cell_id"])
        logger.debug("out: {}".format(out))

        self.assertSequenceEqual(out.vs["id"], self.sym_g.vs["id"])
        self.assertSequenceEqual(out.vs["cell_id"], self.sym_g.vs["cell_id"])
        self.assertSequenceEqual(out.vs["pert_type"], self.sym_g.vs["pert_type"])
        self.assertSequenceEqual(out.es["weight"], self.sym_g.es["weight"])

    def test_asym_gct_to_graph(self):
        logger.debug("self.asym_gct:\n{}".format(self.asym_gct))

        with self.assertRaises(Exception) as e:
            tasseography.asym_gct_to_graph(self.asym_gct, ["pert_time"], ["pert_Time"])
        self.assertIn("field pert_Time not in column metadata", str(e.exception))

        out = tasseography.asym_gct_to_graph(self.asym_gct, ["pert_time"], ["cell_id", "pert_time"])
        logger.debug("out: {}".format(out))

        self.assertSequenceEqual(out.vs["id"], self.asym_g.vs["id"])
        self.assertSequenceEqual(out.vs["type"], self.asym_g.vs["type"])
        self.assertSequenceEqual(out.vs["cell_id"], self.asym_g.vs["cell_id"])
        self.assertSequenceEqual(out.vs["pert_time"], self.asym_g.vs["pert_time"])

        # Use numpy testing to get around NaN
        np.testing.assert_array_equal(out.es["weight"], self.asym_g.es["weight"])

    def test_add_color_attribute_to_vertices(self):
        g_copy = self.sym_g.copy()
        tasseography.add_color_attribute_to_vertices(g_copy, "cell_id")

        # Colors returned don't appear to be stable; OSX v. Linux difference?
        # So just check that it's the same color 3x
        self.assertEqual(g_copy.vs["color"][0], g_copy.vs["color"][1])
        self.assertEqual(g_copy.vs["color"][0], g_copy.vs["color"][2])
        self.assertEqual(g_copy.vs["color"][1], g_copy.vs["color"][2])

        g_copy2 = self.asym_g.copy()
        tasseography.add_color_attribute_to_vertices(g_copy2, "cell_id")

        # Make sure the first 3 colors are the same, that the last 2 colors
        # aren't the same, and that either of the last 2 colors isn't
        # one of the first 3
        self.assertEqual(g_copy2.vs["color"][0], g_copy2.vs["color"][1])
        self.assertEqual(g_copy2.vs["color"][0], g_copy2.vs["color"][2])
        self.assertEqual(g_copy2.vs["color"][1], g_copy2.vs["color"][2])
        self.assertNotEqual(g_copy2.vs["color"][0], g_copy2.vs["color"][3])
        self.assertNotEqual(g_copy2.vs["color"][0], g_copy2.vs["color"][4])
        self.assertNotEqual(g_copy2.vs["color"][3], g_copy2.vs["color"][4])

    def test_add_color_attribute_to_edges(self):
        r = (1., 0., 0., 1.)
        b = (0., 0., 1., 1.)
        k = (0., 0., 0., 1.)

        g_copy = self.sym_g.copy()
        tasseography.add_color_attribute_to_edges(g_copy)
        self.assertSequenceEqual(g_copy.es["color"], [r, r, b])

        g_copy2 = self.asym_g.copy()
        tasseography.add_color_attribute_to_edges(g_copy2)
        self.assertSequenceEqual(g_copy2.es["color"], [r, r, b, b, k, r])

    def test_remove_edges_and_vertices_below_thresh(self):

        out = tasseography.remove_edges_and_vertices_below_thresh(self.asym_g, 0.6, True)
        self.assertSequenceEqual(out.es["weight"], [-0.7, 0.9])
        self.assertSequenceEqual(out.vs["id"], ["B", "C", "o", "k"])

        out2 = tasseography.remove_edges_and_vertices_below_thresh(self.asym_g, 0.6, False)
        self.assertSequenceEqual(out2.es["weight"], [-0.7, 0.9])
        self.assertSequenceEqual(out2.vs["id"], ["A", "B", "C", "o", "k"])

    def test_get_vertex_ids(self):

        out = tasseography.get_vertex_ids(self.sym_g, ["great", "ok"], "pert_type", None)
        self.assertItemsEqual(out, [0, 2])

        out2 = tasseography.get_vertex_ids(self.sym_g, None, "unimportant", None)
        self.assertItemsEqual(out2, [0, 1, 2])

        with self.assertRaises(Exception) as e:
            tasseography.get_vertex_ids(self.sym_g, "a", "unimportant", None)
        self.assertIn("my_query must be a list", str(e.exception))

        out = tasseography.get_vertex_ids(self.asym_g, ["A", "B"], "id", "row")
        self.assertItemsEqual(out, [0, 1])

        out2 = tasseography.get_vertex_ids(self.asym_g, ["A375"], "cell_id", "col")
        self.assertItemsEqual(out2, [4])

    def test_get_vertex_ids_of_neighbors(self):

        out = tasseography.get_vertex_ids_of_neighbors(self.sym_g, [0])
        self.assertItemsEqual(out, set([0, 1, 2]))

        out2 = tasseography.get_vertex_ids_of_neighbors(self.asym_g, [0, 1])
        self.assertItemsEqual(out2, set([0, 1, 3, 4]))

    def test_create_induced_subgraph(self):
        out = self.sym_g.induced_subgraph(set([0, 2]))

        self.assertSequenceEqual(out.vs["id"], ["a", "c"])
        self.assertSequenceEqual(out.es["weight"], [0.6])

        out2 = self.asym_g.induced_subgraph(set([2, 3, 4]))
        self.assertSequenceEqual(out2.vs["id"], ["C", "o", "k"])
        self.assertSequenceEqual(out2.es["weight"], [np.nan, 0.9])

    def test_plot_network(self):
        out_fig_name = "test_plot_network_output.png"

        tasseography.plot_network(self.sym_g, out_fig_name, "id")

        # Make sure plot was produced
        self.assertTrue(os.path.exists(out_fig_name))

        # Remove it
        os.remove(out_fig_name)

    def test_plot_bipartite(self):
        out_fig_name = "test_plot_bipartite_output.png"

        tasseography.plot_network(
            self.asym_g, out_fig_name, "id", layout="bipartite")

        # Make sure plot was produced
        self.assertTrue(os.path.exists(out_fig_name))

        # Remove it
        os.remove(out_fig_name)

    def test_write_graph_to_gml(self):
        out_name = "test_write_graph_to_gml_output.gml"

        tasseography.write_graph_to_gml(self.sym_g, out_name)

        # Make sure text file was produced
        self.assertTrue(os.path.exists(out_name))

        # Remove it
        os.remove(out_name)


if __name__ == '__main__':
    setup_logger.setup(verbose=True)
    unittest.main()
