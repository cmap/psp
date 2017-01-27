"""
network_query_of_gct_igraph.py

Query a connectivity gct and return a network visualization of top connections
using igraph.

"""

import argparse
import logging
import igraph as ig
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_cmap.io.pandasGEXpress.parse as pg

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--in_path", "-i", required=True, help="path to input gct")

    # Optional args
    parser.add_argument("--out_fig_name", "-of", default="network.png",
                        help=("what to name the output figure; file extension " +
                              "will determine file format (e.g. png, svg)"))
    parser.add_argument("--annot_fields", "-a", nargs="*", default=["pert_iname", "moa"],
                        help="fields from metadata to use for annotating nodes")
    parser.add_argument("--threshold", "-t", default=0.8, type=float, help="connectivity threshold")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):

    # Read gct
    gct = pg.parse(args.in_path)

    # Assert column and row metadata fields are the same
    # Assert data_df is square
    assert gct.data_df.shape[0] == gct.data_df.shape[1]

    # Convert to Graph
    g = convert_gct_to_graph(gct, args.annot_fields)

    # Trim edges below the threshold
    trim_graph(g, args.threshold)

    # Write nodes and edges to file

    # # Plot result
    plot_network(g, args.out_fig_name)


def convert_gct_to_graph(gct, annotation_fields):
    """ Convert gct to an igraph.Graph object.

    Args:
        gct (GCToo object)
        annotation_fields (list of strings): metadata to use for annotating nodes

    Returns:
        g (igraph.Graph)

    """

    # Create adjacency matrix
    adj = gct.data_df.values.tolist()

    # Create graph using adjacency matrix
    g = ig.Graph.Weighted_Adjacency(adj, mode=ig.ADJ_DIRECTED,
                                    attr="weight", loops=False)

    # Annotate nodes
    for field in annotation_fields:
        assert field in gct.col_metadata_df.columns.values, (
               "field {} not present in the metadata.".format(field))
        g.vs[field] = gct.col_metadata_df[field]

    return g


def trim_graph(g, thresh):
    """ Remove edges and vertices that are below the threshold. Modify g
    in-place.

    Args:
        g (Graph object)
        thresh (float)

    Returns:
        None

    """

    # Remove edges
    edges_to_delete = [e.index for e in g.es if not e['weight'] > thresh]
    g.delete_edges(edges_to_delete)

    # Remove vertices
    vertices_to_delete = [ind for ind, degree in enumerate(g.degree()) if degree == 0]
    g.delete_vertices(vertices_to_delete)


def plot_network(g, fig_name):

    # Color according to MOA
    unique_moa = set(g.vs["moa"])
    num_moa = len(unique_moa)
    palette = ig.drawing.colors.ClusterColoringPalette(num_moa)
    moa_color_dict = dict(zip(unique_moa, palette))
    moa_colors = [moa_color_dict[moa] for moa in g.vs["moa"]]

    ig.plot(g, fig_name,
            vertex_size=10,
            vertex_label=g.vs['pert_iname'],
            vertex_color=moa_colors,
            vertex_label_dist=2,
            edge_arrow_size=0.3,
            edge_curved=False)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)