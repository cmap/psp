"""
tasseography.py

Collection of functions that produce networks from gct files. Includes
functionality for both symmetric and asymmetric gcts.

"""

import logging
import sys
import argparse
import numpy as np
import igraph as ig

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
    parser.add_argument("--my_query", "-q", nargs="*", help="query values")
    parser.add_argument("--annot_fields", "-a", nargs="*", default=["pert_iname", "moa"],
                        help="fields from metadata to use for annotating vertices")
    parser.add_argument("--threshold", "-t", default=0.8, type=float,
                        help="connectivity threshold")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):

    # Parse gct
    gct = pg.parse(args.in_path)

    # Figure out whether or not the gct is symmetric
    if gct.row_metadata_df.equals(gct.col_metadata_df):
        logger.info("Row metadata equals column metadata. Assuming symmetric GCT.")
        g = sym_gct_to_graph(gct, annot_fields)

    else:
        logger.info("Row metadata does not equal column metadata. Assuming asymmetric GCT.")
        g = asym_gct_to_graph(gct, row_annot_fields, col_annot_fields)

    # # Remove edges below thresh
    # remove_edges_below_thresh(g, args.threshold)
    #
    # # Get ids of vertices to keep
    # get_ids_of_vertices_sym()




def main_sym(gct, out_fig_name, vertex_annot_fields, my_query,
             my_query_annot_field, threshold, vertex_label_field,
             vertex_color_field, layout):

    # Convert gct to Graph object
    g = sym_gct_to_graph(gct, vertex_annot_fields)

    # TODO(LL): make color dict

    # Remove edges (and optionally vertices) from subgraph below threshold
    subgraph = remove_edges_and_vertices_below_thresh(g, threshold)

    # Get vertex ids for my_query
    vertex_ids_of_queries = get_vertex_ids_sym(subgraph, my_query, my_query_annot_field)

    # Add in vertex ids for all first-order neighbors
    vertex_ids_of_queries_and_neighbors = get_vertex_ids_of_neighbors(
        subgraph, vertex_ids_of_queries)
    logger.info("Graph has {} vertices.".format(len(vertex_ids_of_queries_and_neighbors)))

    # Keep only vertices and edges related to queries and their neighbors
    out_graph = subgraph.induced_subgraph(vertex_ids_of_queries_and_neighbors)

    # TODO(LL): export GML

    # Plot network if out_fig_name provided
    if out_fig_name:
        plot_network(out_graph, out_fig_name,
                     vertex_label_field=vertex_label_field,
                     vertex_color_field=vertex_color_field,
                     layout=layout)


def main_asym(gct, out_fig_name, row_annot_fields, col_annot_fields,
              my_query_in_rows, my_query_in_cols,
              my_query_annot_field_in_rows, my_query_annot_field_in_cols,
              threshold, vertex_label_field, vertex_color_field):

    # Convert gct to Graph object
    g = asym_gct_to_graph(gct, row_annot_fields, col_annot_fields)

    # TODO(LL): make color dict

    # Remove edges (and optionally vertices) from subgraph below threshold
    subgraph = remove_edges_and_vertices_below_thresh(g, threshold)

    # Get vertex ids for my_query_in_rows and my_query_in_cols
    vertex_ids_of_queries = get_vertex_ids_asym(
        subgraph, my_query_in_rows, my_query_in_cols,
        my_query_annot_field_in_rows, my_query_annot_field_in_cols)

    # TODO(LL): the logic needs to reworked; right now, there is no way to ask
    # "show me all connections above thresh for these elements in the columns
    # to any elements in the rows"

    # Add in vertex ids for all first-order neighbors
    vertex_ids_of_queries_and_neighbors = get_vertex_ids_of_neighbors(
        subgraph, vertex_ids_of_queries)
    logger.info("Graph has {} vertices.".format(len(vertex_ids_of_queries_and_neighbors)))

    # Keep only vertices and edges related to queries and their neighbors
    out_graph = subgraph.induced_subgraph(vertex_ids_of_queries_and_neighbors)

    # TODO(LL): export GML

    # Plot bipartite graph if out_fig_name provided
    if out_fig_name:
        plot_network(out_graph, out_fig_name,
                     vertex_label_field=vertex_label_field,
                     vertex_color_field=vertex_color_field,
                     layout="bipartite")


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


def get_vertex_ids_sym(g, my_query, my_query_annot_field):
    """ Extract vertices with the values in my_query in my_query_annot_field.

    Args:
        @param g:
        @type g: igraph.Graph object
        @param my_query: vertices to keep
        @type my_query: list of strings or None
        @param my_query_annot_field: vertex attribute field in which to look for my_query
        @param my_query_annot_field: string

    Returns:
        vertex_ids (list of integers): vertex ids corresponding to my_query

    """
    # If no query provided, return all vertex ids
    if my_query is None:
        vertex_ids = [v.index for v in g.vs()]

    elif type(my_query) is list:

        # Find all vertices matching the values in my_query
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


def get_vertex_ids_asym(g, my_query_in_rows, my_query_in_cols,
                        my_query_annot_field_in_rows, my_query_annot_field_in_cols):
    """ Extract vertices from the bipartite graph g. Row and column vertices are
    specified separately.

    Args:
        @param g
        @type g: igraph.Graph object
        @param my_query_in_rows: vertices to keep in rows
        @type my_query_in_rows: list of strings or None
        @param my_query_in_cols: vertices to keep in columns
        @type my_query_in_cols: list of strings or None
        @param my_query_annot_field_in_rows: vertex attribute field in which to look for my_query_in_rows
        @param my_query_annot_field_in_rows: string
        @param my_query_annot_field_in_cols: vertex attribute field in which to look for my_query_in_cols
        @param my_query_annot_field_in_cols: string

    Returns:
        vertex_ids (list of integers): vertex ids corresponding to my_query_in_cols and my_query_in_rows

    """
    # Initialize vertex_ids
    vertex_ids = []

    ### First, get row vertices

    # If no query provided, return all vertex ids; v["type"] == False for rows
    if my_query_in_rows is None:
        vertex_ids += [v.index for v in g.vs() if not v["type"]]

    elif type(my_query_in_rows) is list:

        # Find all row vertices matching the values in my_query_in_rows
        for row_q in my_query_in_rows:
            these_ids = [v.index for v in g.vs() if not v["type"] and v[my_query_annot_field_in_rows] == row_q]
            if len(these_ids) < 1:
                msg = "No vertices were found for query {} in the row vertex attribute {}.".format(row_q, my_query_annot_field_in_rows)
                logger.info(msg)
            vertex_ids += these_ids

    else:
        msg = "my_query_in_rows must be a list. my_query_in_rows: {}".format(my_query_in_rows)
        logger.error(msg)
        raise Exception(msg)

    ### Now, get column vertices

    # If no query provided, return all vertex ids; v["type"] == True for columns
    if my_query_in_cols is None:
        vertex_ids += [v.index for v in g.vs() if v["type"]]

    elif type(my_query_in_cols) is list:

        # Find all column vertices matching the values in my_query_in_cols
        for col_q in my_query_in_cols:
            these_ids = [v.index for v in g.vs() if v["type"] and v[my_query_annot_field_in_cols] == col_q]
            if len(these_ids) < 1:
                msg = "No vertices were found for query {} in the column vertex attribute {}.".format(col_q, my_query_annot_field_in_cols)
                logger.info(msg)
            vertex_ids += these_ids

    else:
        msg = "my_query_in_cols must be a list. my_query_in_cols: {}".format(my_query_in_cols)
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


def plot_network(g, out_fig_name, vertex_label_field, vertex_color_field,layout="fr"):
    """ Plot network.

    Args:
        @param g
        @type g: iGraph.Graph object
        @param out_fig_name: if None, image won't be saved to file
        @type out_fig_name: string or None
        @param vertex_label_field: vertex attribute field used for labels
        @type vertex_label_field: string
        @param vertex_color_field: vertex attribute field used for colors
        @type vertex_color_field: string
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

    # Check if vertex_color_field was provided
    if vertex_color_field is None:
        colors = None
    else:
        colors = g.vs[vertex_color_field]

    # Do it
    ig.plot(g, out_fig_name,
            layout=layout,
            vertex_size=10,
            vertex_label=labels,
            vertex_color=colors,
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
    g.write_gml(out_gml_name)
    logger.info("Graph written to {}.".format(out_gml_name))


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)