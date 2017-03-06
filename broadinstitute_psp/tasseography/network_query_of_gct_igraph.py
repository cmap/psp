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
    parser.add_argument("--out_node_name", "-on", default="nodes.tsv",
                        help="what to name output nodes file")
    parser.add_argument("--out_edge_name", "-oe", default="edges.tsv",
                        help="what to name output edges file")
    parser.add_argument("--annot_fields", "-a", nargs="*",
                        default=["pert_iname", "moa"],
                        help="fields from metadata to use for annotating nodes")
    parser.add_argument("--node_color_field", "-nc", default="moa",
                        help="field from metadata to use for coloring nodes")
    parser.add_argument("--threshold", "-t", default=0.8, type=float, help="connectivity threshold")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


# TODO(lev): make a (bipartite) version of this script for assymetric matrices


def main(args):

    # Read gct
    gct = pg.parse(args.in_path)

    # Assert column and row metadata fields are the same
    # Assert data_df is square
    assert gct.data_df.shape[0] == gct.data_df.shape[1]

    # Make dictionary of node colors for plotting
    color_dict = make_color_dict(gct.col_metadata_df, args.node_color_field)

    # Convert to Graph
    g = convert_gct_to_graph(gct, args.annot_fields)

    # Trim edges below the threshold
    trim_graph(g, args.threshold)

    # Write nodes and edges to file
    write_edge_tsv(g, args.out_edge_name)
    write_node_tsv(g, args.out_node_name)

    # Plot result
    plot_network(g, args.out_fig_name, color_dict, args.node_color_field)


def convert_gct_to_graph(gct, annotation_fields):
    """ Convert gct to an igraph.Graph object. Id is annotated as the id
    attribute.

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

    # Annotate nodes using ids
    g.vs["id"] = gct.col_metadata_df.index.values

    # Annotate nodes with other metadata
    for field in annotation_fields:
        assert field in gct.col_metadata_df.columns.values, (
               "field {} not present in the metadata.".format(field))
        g.vs[field] = gct.col_metadata_df[field]

    return g


def make_color_dict(meta_df, color_field):
    """ Keys are unique entries in df[color_field], values are RGB color tuples.

    Args:
        meta_df (pandas df)
        color_field (string): metadata field to use for node colors

    Returns:
        color_dict

    """
    unique_vals = meta_df[color_field].unique()
    num_unique_vals = len(unique_vals)

    palette = ig.drawing.colors.ClusterColoringPalette(num_unique_vals)
    color_dict = dict(zip(unique_vals, palette))

    return color_dict


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
    logger.info("Threshold: abs(e['weight']) > {thresh}]".format(thresh=thresh))
    edges_to_delete = [e.index for e in g.es if not abs(e['weight']) > thresh]
    g.delete_edges(edges_to_delete)

    # Remove vertices
    vertices_to_delete = [ind for ind, degree in enumerate(g.degree()) if degree == 0]
    g.delete_vertices(vertices_to_delete)


def write_edge_tsv(g, out_edge_name):
    """ Write edges to a tsv.

    Args:
        g (igraph.Graph)
        out_edge_name (file path)

    Returns:
        None

    """
    # Make dataframe
    source_target_weight_tuples = [(g.vs[e.source]["id"], g.vs[e.target]["id"],
                                    e["weight"]) for e in g.es]
    edge_df = pd.DataFrame(source_target_weight_tuples,
                           columns=["query", "target", "value"])

    # Write to file
    edge_df.to_csv(out_edge_name, sep="\t", index=None)


def write_node_tsv(g, out_node_name):
    """ Write nodes (or vertices) to a tsv.

    Args:
        g (igraph.Graph)
        meta_df (pandas df)
        annot_fields (list of strings)
        out_node_name (file path)

    Returns:
        None

    """
    # Make dataframe
    node_attributes = g.vs.attribute_names()
    node_df = pd.DataFrame.from_dict(
        dict([(attr, g.vs[attr]) for attr in node_attributes]))

    # Write to file
    node_df.to_csv(out_node_name, sep="\t", index=None)


def plot_network(g, fig_name, color_dict, color_field, layout="fr"):
    margin = 80

    node_colors = [color_dict[val] for val in g.vs[color_field]]
    ig.plot(g, fig_name,
            layout=layout,
            vertex_size=10,
            vertex_label=g.vs['pert_iname'],
            vertex_color=node_colors,
            vertex_label_dist=2,
            edge_arrow_size=0.3,
            edge_curved=False,
            margin=(margin, margin, margin, margin))

    logger.info("Network saved to {}".format(fig_name))


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)