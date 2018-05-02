"""
tasseography.py

Collection of functions that produce networks from gct files. Networks are
trimmed according to the provided threshold; only edges with an absolute value
above the threshold are returned. If only a subgraph is desired,
my_query can be used to specify one or more vertices. Only those vertices and
their first-order neighbors will be returned in the figure and gml file.

This module includes functionality for both symmetric and asymmetric gcts.
In the case of an asymmetric gct, the output figure will be a bipartite graph.

There are lots of arguments; the most important ones are input_gct_path,
threshold, and my_query.

"""

import logging
import sys
import argparse
import numpy as np
import igraph as ig
import warnings
import matplotlib.pyplot as plt

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.parse as parse

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Default layout is Fruchterman-Reingold
LAYOUT = "fruchterman_reingold"
COLORMAP = plt.cm.gnuplot

def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--input_gct_path", "-i", required=True,
                        help="path to input gct")

    # Optional args
    parser.add_argument("--out_fig_name", "-of", default=None,
                        help=("what to name the output figure; if not provided, " +
                              "no figure is produced; extension determines " +
                              "file type"))
    parser.add_argument("--out_gml_name", "-og", default=None,
                        help=("what to name the output gml file; if not provided, " +
                              "no gml file is produced; GML is a common text file " +
                              "format for graphs"))

    meg = parser.add_mutually_exclusive_group()
    meg.add_argument("--threshold", "-t", type=float,
                        help=("only edges whose ABSOLUTE value is above the " +
                              "threshold will be returned"))
    meg.add_argument("--percentile", "-p", type=float,
                        help=("only edges whose ABSOLUTE value is above this " +
                              "percentile threshold will be returned; " +
                              "provide either percentile or threshold, " +
                              "but not both"))

    parser.add_argument("--my_query", "-q", nargs="*", default=None,
                        help=("only subgraphs including my_query will be " +
                              "included; set to None to see all " +
                              "connections above threshold"))
    parser.add_argument("--query_field", "-qf", type=str,
                        default="pert_iname",
                        help=("field from metadata in which to look for " +
                              "my_query; use 'id' to look in the indices"))
    parser.add_argument("--query_in_row_or_col", "-rc", default=None,
                        choices=["row", "col", None],
                        help=("indicates whether to look for my_query in " +
                              "query_field in the row or column metadata; " +
                              "set to None for symmetric GCTs;"))

    parser.add_argument("--row_annot_fields", "-ra", nargs="*",
                        default=["pert_iname", "moa", "cell_id"],
                        help=("fields from row metadata to use for " +
                              "annotating row vertices"))
    parser.add_argument("--col_annot_fields", "-ca", nargs="*",
                        default=["pert_iname", "moa", "cell_id"],
                        help=("fields from column metadata to use for " +
                              "annotating column vertices; should be the " +
                              "as row_annot_fields if GCT is symmetric"))
    parser.add_argument("--vertex_label_field", "-vl", default=None,
                        help=("metadata field to use for labeling " +
                              "vertices in the figure"))
    parser.add_argument("--vertex_color_field", "-vc", default=None,
                        help=("metadata field to use for coloring " +
                              "vertices in the figure"))
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):

    # Parse gct
    gct = parse.parse(args.input_gct_path)

    # TODO(LL): better integrate main_sym and main_asym

    # Figure out whether or not the gct is symmetric
    if gct.row_metadata_df.equals(gct.col_metadata_df):
        logger.info(("Row metadata equals column metadata. " +
                     "Assuming symmetric GCT."))

        assert args.row_annot_fields == args.col_annot_fields, (
            ("row_annot_fields should be the same as col_annot_fields if the " +
             "GCT is symmetric. args.row_annot_fields: {}, " +
             "args.col_annot_fields: {}").format(args.row_annot_fields,
                                                 args.col_annot_fields))

        assert args.query_in_row_or_col is None, (
            ("query_in_row_or_col should be None for symmetric GCTs. " +
             "args.query_in_row_or_col: {}").format(args.query_in_row_or_col))

        # Main method for symmetric gcts
        main_sym(gct, args.out_fig_name, args.out_gml_name,
                 args.row_annot_fields, args.my_query, args.query_field,
                 args.threshold, args.percentile, args.vertex_label_field,
                 args.vertex_color_field, layout=LAYOUT)

    else:
        logger.info(("Row metadata does not equal column metadata. " +
                     "Assuming asymmetric GCT."))

        assert args.query_in_row_or_col != "both", (
            ("query_in_row_or_col must not be 'both' if the matrix is " +
             "asymmetric. args.query_in_row_or_col: {}").format(
                args.query_in_row_or_col))

        # Main method for asymmetric gcts
        main_asym(gct, args.out_fig_name, args.out_gml_name,
                  args.row_annot_fields, args.col_annot_fields,
                  args.my_query, args.query_field, args.query_in_row_or_col,
                  args.threshold, args.percentile, args.vertex_label_field,
                  args.vertex_color_field)


def convert_percentile_to_thresh(g, percentile):
    """ Figure out what value in g corresponds to the given percentile.

    Args:
        @param g:
        @type g: igraph.Graph
        @param percentile: between 0 and 100
        @type percentile: float

    Returns:
        thresh

    """
    assert percentile < 100 and percentile > 0, (
        "percentile must be between 0 and 100. percentile: {}".format(percentile))

    thresh = np.nanpercentile(np.abs(g.es["weight"]), percentile)

    return thresh


def main_sym(gct, out_fig_name, out_gml_name, vertex_annot_fields, my_query,
             my_query_annot_field, threshold, percentile, vertex_label_field,
             vertex_color_field, layout):

    # Convert gct to Graph object
    g = sym_gct_to_graph(gct, vertex_annot_fields)

    # Calculate threshold from percentile if percentile provided
    if percentile is not None:
        thresh = convert_percentile_to_thresh(g, percentile)
    else:
        thresh = threshold

    # Add 'color' attribute to nodes based on entries in vertex_color_field
    if vertex_color_field is not None:
        add_color_attribute_to_vertices(g, vertex_color_field)

    # Add 'color' attribute to edges based on whether 'weight' is + or -
    add_color_attribute_to_edges(g)

    # Remove edges from subgraph below threshold
    subgraph = remove_edges_and_vertices_below_thresh(g, thresh)

    # Get vertex ids for my_query
    vertex_ids_of_queries = get_vertex_ids(subgraph, my_query,
                                           my_query_annot_field, None)

    # Add in vertex ids for all first-order neighbors
    vertex_ids_of_queries_and_neighbors = get_vertex_ids_of_neighbors(
        subgraph, vertex_ids_of_queries)

    # Keep only vertices and edges related to queries and their neighbors
    out_graph = subgraph.induced_subgraph(vertex_ids_of_queries_and_neighbors)

    logger.info("Graph has {} vertices.".format(out_graph.vcount()))
    logger.info("Graph has {} edges.".format(out_graph.ecount()))

    # Write graph to .gml file if out_gml_name provided
    if out_gml_name:
        write_graph_to_gml(out_graph, out_gml_name)

    # Plot network if out_fig_name provided
    if out_fig_name:
        plot_network(out_graph, out_fig_name,
                     vertex_label_field=vertex_label_field,
                     layout=layout)


def main_asym(gct, out_fig_name, out_gml_name, row_annot_fields, col_annot_fields,
              my_query, my_query_annot_field, query_in_row_or_col,
              threshold, percentile, vertex_label_field, vertex_color_field):

    # Convert gct to Graph object
    g = asym_gct_to_graph(gct, row_annot_fields, col_annot_fields)

    # Calculate threshold from percentile if percentile provided
    if percentile is not None:
        thresh = convert_percentile_to_thresh(g, percentile)
    else:
        thresh = threshold

    # Add 'color' field to use for coloring vertices
    if vertex_color_field is not None:
        add_color_attribute_to_vertices(g, vertex_color_field)

    # Add 'color' attribute to edges
    add_color_attribute_to_edges(g)

    # Remove edges (and optionally vertices) from subgraph below threshold
    subgraph = remove_edges_and_vertices_below_thresh(g, thresh)

    # Get vertex ids for my_query_in_rows and my_query_in_cols
    vertex_ids_of_queries = get_vertex_ids(
        subgraph, my_query, my_query_annot_field, query_in_row_or_col)

    # Add in vertex ids for all first-order neighbors
    vertex_ids_of_queries_and_neighbors = get_vertex_ids_of_neighbors(
        subgraph, vertex_ids_of_queries)
    logger.info("Graph has {} vertices.".format(len(vertex_ids_of_queries_and_neighbors)))

    # Keep only vertices and edges related to queries and their neighbors
    out_graph = subgraph.induced_subgraph(vertex_ids_of_queries_and_neighbors)

    # Write graph to .gml file if out_gml_name provided
    if out_gml_name:
        write_graph_to_gml(out_graph, out_gml_name)

    # Plot bipartite graph if out_fig_name provided
    if out_fig_name:

        # Rotating the layout helps a lot
        bipartite_layout = out_graph.layout_bipartite()
        bipartite_layout.rotate(-90)

        plot_network(out_graph, out_fig_name,
                     vertex_label_field=vertex_label_field,
                     layout=bipartite_layout)


def sym_gct_to_graph(gct, annot_fields):
    """ Convert symmetric gct to an igraph.Graph object. Row indices are used to
    populate the "id" attribute of the vertices of the graph. The maximum value
    of A(i,j), A(j,i) is used.

    Args:
        @param gct:
        @type gct: GCToo object
        @param annot_fields: metadata to use for annotating vertices
        @type annot_fields: list of strings

    Returns:
        g (igraph.Graph)

    """
    assert gct.row_metadata_df.equals(gct.col_metadata_df), (
        "Row metadata must be the same as the column metadata.")

    # Create adjacency matrix
    adj = gct.data_df.values.tolist()

    # Create graph using adjacency matrix
    g = ig.Graph.Weighted_Adjacency(adj, mode=ig.ADJ_MAX, attr="weight", loops=False)

    # Annotate vertices using ids
    g.vs["id"] = gct.row_metadata_df.index.values

    # Annotate vertices
    for field in annot_fields:
        assert field in gct.row_metadata_df.columns.values, (
               "field {} not present in the metadata.".format(field))
        g.vs[field] = gct.row_metadata_df[field]

    return g


def asym_gct_to_graph(gct, row_annot_fields, col_annot_fields):
    """ Convert asymmetric gct to an igraph.Graph object. Annotations in the
    rows are kept separate from the annotations in the columns.

    Args:
        @param gct:
        @type gct: GCToo object
        @param row_annot_fields: metadata to use for annotating row vertices
        @type row_annot_fields: list of strings
        @param col_annot_fields: metadata to use for annotating col vertices
        @type col_annot_fields: list of strings

    Returns:
        g (igraph.Graph)

    """
    # Initialize graph
    g = ig.Graph.Full_Bipartite(gct.data_df.shape[0], gct.data_df.shape[1])

    # Assign weights to each edge
    g.es["weight"] = gct.data_df.values.flatten()

    # Assign id to each vertex (row vertices are first, then column vertices)
    g.vs["id"] = np.concatenate([gct.data_df.index.values, gct.data_df.columns.values])

    # Add annotations
    for v in g.vs():

        # v["type"] = True implies column vertex
        if v["type"]:

            # Get column metadata
            for annot_field in col_annot_fields:
                assert annot_field in gct.col_metadata_df.columns.values, (
                    ("field {} not in column metadata. gct.col_metadata_df." +
                    "columns.values: {}").format(
                        annot_field, gct.col_metadata_df.columns.values))
                v[annot_field] = gct.col_metadata_df.loc[v["id"], annot_field]

        # v["type"] = False implies row vertex
        else:

            # Get row metadata
            for annot_field in row_annot_fields:
                assert annot_field in gct.row_metadata_df.columns.values, (
                    ("field {} not in row metadata. gct.row_metadata_df." +
                     "columns.values: {}").format(
                        annot_field, gct.row_metadata_df.columns.values))
                v[annot_field] = gct.row_metadata_df.loc[v["id"], annot_field]

    return g


def add_color_attribute_to_vertices(g, vertex_color_field):
    """ Add a vertex attribute called 'color' that corresponds to an existing
     vertex attribute vertex_color_field.

    Args:
        @param g:
        @type g: iGraph.Graph object
        @param vertex_color_field: name of existing vertex attribute
        @type vertex_color_field: string

    Returns:
        None

    """
    # Figure out how many unique elements we have in vertex_color_field
    unique_elements = set(g.vs[vertex_color_field])

    # Create a dictionary in which each unique element is assigned a color
    color_dict = dict(zip(
        unique_elements, [COLORMAP(jj) for jj in np.linspace(0, 0.9, len(unique_elements))]))

    # Populate the 'color' vertex attribute
    list_of_colors = [color_dict[v[vertex_color_field]] for v in g.vs()]
    g.vs["color"] = list_of_colors


def add_color_attribute_to_edges(g):
    """ Add an edge attribute called 'color' indicating whether the edge has
    positive or negative weight. -1 means negative, 1 means positive, 0 means
    the weight was NaN.

    Args:
        @param g
        @type g: iGraph.Graph object

    Returns:

    """
    POSITIVE_EDGE_COLOR = (1., 0., 0., 1.)
    NEGATIVE_EDGE_COLOR = (0., 0., 1., 1.)
    OTHER_EDGE_COLOR = (0., 0., 0., 1.)

    g.es["color"] = [
        NEGATIVE_EDGE_COLOR if e["weight"] < 0
        else POSITIVE_EDGE_COLOR if e["weight"] > 0
        else OTHER_EDGE_COLOR for e in g.es()]


def remove_edges_and_vertices_below_thresh(g, thresh, delete_vertices=True):
    """ Remove edges and vertices that are below the threshold.

    Args:
        @param g:
        @type g: igraph.Graph object
        @param thresh:
        @type thresh: float
        @param delete_vertices: whether to delete vertices without edges
        @type delete_vertices: bool

    Returns:
        subgraph (igraph.Graph object)

    """
    # Remove edges
    logger.info("Threshold: abs(e['weight']) > {thresh}]".format(thresh=thresh))
    edges_to_keep = [e.index for e in g.es if abs(e['weight']) > thresh]
    subgraph = g.subgraph_edges(edges_to_keep, delete_vertices)

    return subgraph


def get_vertex_ids(g, my_query, my_query_annot_field, row_or_col):
    """ Extract vertices with the values in my_query in my_query_annot_field.
    If row_or_col is "row", my_query will only be searched for in row vertices.
    If row_or_col is "col", my_query will only be searched for in column
    vertices. If row_or_col is None, my_query will be searched for in all
    vertices.

    Args:
        @param g:
        @type g: igraph.Graph object
        @param my_query: vertices to keep
        @type my_query: list of strings or None
        @param my_query_annot_field: vertex attribute field in which to look for my_query
        @param my_query_annot_field: string
        @param row_or_col: indicates which vertices to look in
        @type row_or_col: string or id

    Returns:
        vertex_ids (list of integers): vertex ids corresponding to my_query

    """
    assert row_or_col in ["row", "col", None], (
        "row_or_col must be 'row', 'col', or None. row_or_col: {}".format(row_or_col))

    # If no query provided, return all vertex ids
    if my_query is None:
        if row_or_col == "row":
            vertex_ids = [v.index for v in g.vs() if not v["type"]]
        elif row_or_col == "col":
            vertex_ids = [v.index for v in g.vs() if v["type"]]
        else:
            vertex_ids = [v.index for v in g.vs()]

    elif type(my_query) is list:

        if row_or_col == "row":

            # Find all ROW vertices matching the values in my_query
            vertex_ids = []
            for q in my_query:
                these_ids = [v.index for v in g.vs() if v[my_query_annot_field] == q and not v["type"]]
                if len(these_ids) < 1:
                    msg = "No vertices were found for query {} in row vertex attribute {}.".format(q, my_query_annot_field)
                    logger.info(msg)
                vertex_ids += these_ids

        elif row_or_col == "col":

            # Find all COLUMN vertices matching the values in my_query
            vertex_ids = []
            for q in my_query:
                these_ids = [v.index for v in g.vs() if v[my_query_annot_field] == q and v["type"]]
                if len(these_ids) < 1:
                    msg = "No vertices were found for query {} in column vertex attribute {}.".format(q, my_query_annot_field)
                    logger.info(msg)
                vertex_ids += these_ids

        else:

            # Find ALL vertices matching the values in my_query
            vertex_ids = []
            for q in my_query:
                these_ids = [v.index for v in g.vs() if v[my_query_annot_field] == q]
                if len(these_ids) < 1:
                    msg = "No vertices were found for query {} in vertex attribute {}.".format(q, my_query_annot_field)
                    logger.info(msg)
                vertex_ids += these_ids

    else:
        msg = "my_query must be a list. my_query: {}".format(my_query)
        logger.error(msg)
        raise Exception(msg)

    return vertex_ids


def get_vertex_ids_of_neighbors(g, vertex_ids_of_queries):
    """ Return first-order neighbors for the provided vertex_ids. Result
    will include the input vertex_ids in addition to the neighbors.

    Args:
        @param g
        @type g: igraph.Graph object
        @param vertex_ids_of_queries:
        @type vertex_ids_of_queries: list of integers

    Returns:
        vertex_ids_of_neighbors (set of integers)

    """
    # Return neighborhood around each query
    list_of_lists_of_neighbors = g.neighborhood(vertex_ids_of_queries)

    # Convert the list of lists to a set
    vertex_ids_of_neighbors = set([vertex_id for sublist in list_of_lists_of_neighbors for vertex_id in sublist])

    return vertex_ids_of_neighbors


def plot_network(g, out_fig_name, vertex_label_field, layout="fr"):
    """ Plot network.

    Args:
        @param g
        @type g: iGraph.Graph object
        @param out_fig_name: if None, image won't be saved to file
        @type out_fig_name: string or None
        @param vertex_label_field: vertex attribute field used for labels
        @type vertex_label_field: string
        @param layout: layout to use for plotting; default = "fr", aka
            Fruchterman Reingold; some other options include "kk" for
            Kamada Kawai, "circle" for circular, and "bipartite"
        @type layout: string

    Returns:
        None

    """
    margin = 80

    # Check if vertex_label_field was provided
    if vertex_label_field is None:
        labels = None
    else:
        labels = g.vs[vertex_label_field]

    # # Check if vertex_color_field was provided
    # if vertex_color_field is None:
    #     colors = None
    # else:
    #     colors = g.vs[vertex_color_field]

    # Do it
    ig.plot(g, out_fig_name,
            layout=layout,
            vertex_size=10,
            vertex_label=labels,
            vertex_label_dist=2,
            edge_arrow_size=0.3,
            edge_curved=False,
            margin=(margin, margin, margin, margin))

    if out_fig_name is not None:
        logger.info("Network saved to {}".format(out_fig_name))


def write_graph_to_gml(g, out_gml_name):
    """ Write graph to .gml (Graph Markup Language) file.

    Args:
        @param g:
        @type g: iGraph.Graph object
        @param out_gml_name: name of output .gml file
        @type out_gml_name: string

    Returns:
        None

    """
    # Ignore warnings thrown by write_gml
    warnings.simplefilter("ignore")

    g.write_gml(out_gml_name)
    logger.info("Graph written to {}.".format(out_gml_name))


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)
