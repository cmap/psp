"""
<NAME_OF_SCRIPT>.py

Describe what this script does.

"""

import logging
import sys
import argparse
import numpy as np
import pandas as pd
import igraph as ig

import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_cmap.io.pandasGEXpress.GCToo as GCToo
import broadinstitute_cmap.io.pandasGEXpress.parse as pg
import broadinstitute_cmap.io.pandasGEXpress.write_gct as wg

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--arg1", "-a1", required=True, help="argument 1")

    # Optional args
    parser.add_argument("--arg2", "-a2", help="argument 2")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):

    # Parse gct
    gct = pg.parse(args.input_gct_path)

    # Figure out whether or not the gct is symmetric
    if gct.row_metadata_df.equals(gct.col_metadata_df):
        g = sym_gct_to_graph(gct, annot_fields)

    else:
        g = asym_gct_to_graph(gct, row_annot_fields, col_annot_fields)






def sym_gct_to_graph(gct, annot_fields):
    """ Convert symmetric gct to an igraph.Graph object. Row indices are used to
    populate the "id" attribute of the nodes of the graph. The maximum value
    of A(i,j), A(j,i) is used.

    Args:
        @param: gct (GCToo object)
        @param: annot_fields (list of strings): metadata to use for annotating nodes

    Returns:
        g (igraph.Graph)

    """
    assert gct.row_metadata_df.equals(gct.col_metadata_df), (
        "Row metadata must be the same as the column metadata.")

    # Create adjacency matrix
    adj = gct.data_df.values.tolist()

    # Create graph using adjacency matrix
    g = ig.Graph.Weighted_Adjacency(adj, mode=ig.ADJ_MAX, attr="weight", loops=False)

    # Annotate nodes using ids
    g.vs["id"] = gct.row_metadata_df.index.values

    # Annotate nodes
    for field in annot_fields:
        assert field in gct.row_metadata_df.columns.values, (
               "field {} not present in the metadata.".format(field))
        g.vs[field] = gct.row_metadata_df[field]

    return g


def asym_gct_to_graph(gct, row_annot_fields, col_annot_fields):
    """ Convert asymmetric gct to an igraph.Graph object. Annotations in the
    rows are kept separate from the annotations in the columns.

    Args:
        @param: gct (GCToo object)
        @param: row_annot_fields (list of strings): metadata to use for annotating row nodes
        @param: col_annot_fields (list of strings): metadata to use for annotating col nodes

    Returns:
        g (igraph.Graph)

    """
    # Initialize graph
    g = ig.Graph.Full_Bipartite(gct.data_df.shape[0], gct.data_df.shape[1])

    # Assign weights to each edge
    g.es["weight"] = gct.data_df.values.flatten()

    # Assign id to each node (row nodes are first, then column nodes)
    g.vs["id"] = np.concatenate([gct.data_df.index.values, gct.data_df.columns.values])

    # Add annotations
    for v in g.vs():

        # v["type"] = True implies column node
        if v["type"]:

            # Get column metadata
            for annot_field in col_annot_fields:
                assert annot_field in gct.col_metadata_df.columns.values, (
                    ("field {} not in column metadata. gct.col_metadata_df." +
                    "columns.values: {}").format(
                        annot_field, gct.col_metadata_df.columns.values))
                v[annot_field] = gct.col_metadata_df.loc[v["id"], annot_field]

        # v["type"] = False implies row node
        else:

            # Get row metadata
            for annot_field in row_annot_fields:
                assert annot_field in gct.row_metadata_df.columns.values, (
                    ("field {} not in row metadata. gct.row_metadata_df." +
                    "columns.values: {}").format(
                        annot_field, gct.row_metadata_df.columns.values))
                v[annot_field] = gct.row_metadata_df.loc[v["id"], annot_field]

    return g

def remove_edges_below_thresh(g, thresh):
    pass

def extract_subgraph(g, node_field, node_value):
    pass

def plot_network(g, layout, node_color_field, edge_color_field):
    pass

def plot_bipartite(g):
    pass

def graph_to_text_files(g, out_node_name, out_edge_name):
    pass



if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)