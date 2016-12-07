"""
network_query_of_gct.py

Produce 2 text files that can be used by Cytoscape and other network
visualizations: a nodes tsv file and an edges tsv file.

my_query can be a single query, a list of queries, or None; if None, all
connections with an absolute value greater or equal to the threshold are
returned.

"""

import logging
import sys
import argparse
import pandas as pd
import numpy as np
import itertools

import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_cmap.io.GCToo.parse_gctoo as pg

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

# Turn off SettingWithCopyWarning
pd.options.mode.chained_assignment = None

logger = logging.getLogger(setup_logger.LOGGER_NAME)

# TODO(lev): make sure annotations are being pulled from row metadata if query
# and from col metadata if target; right now, just pulled only from col metadata


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--in_path", "-i", required=True, help="path to input gct")

    # Optional args
    parser.add_argument("--my_query", "-q", nargs="*", help="pert_iname(s) of query")
    parser.add_argument("--annot_fields", "-a", nargs="*", default=["pert_iname", "moa"],
                        help="fields from metadata to use for annotating nodes")
    parser.add_argument("--out_node_name", "-on", default="nodes.tsv", help="name of output node file")
    parser.add_argument("--out_edge_name", "-oe", default="edges.tsv", help="name of output edge file")
    parser.add_argument("--threshold", "-t", default=0.8, type=float, help="connectivity threshold")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):

    # Read gct
    gct = pg.parse(args.in_path)

    assert all(x in gct.col_metadata_df.columns for x in args.annot_fields), (
        "args.annot_fields {} must be in gct.col_metadata_df.columns: {}".format(
            args.annot_fields, gct.col_metadata_df.columns.values))

    # Melt gct
    melted_df_w_self_conns = melt_gct(gct)

    # Remove self connections
    melted_df = remove_self_cxns(melted_df_w_self_conns)

    # Create edge and node dataframes -- main method
    node_df, edge_df = network_query_of_df(melted_df, gct.col_metadata_df,
                                           args.my_query, args.threshold, args.annot_fields)

    # Write to tsv files
    edge_df.to_csv(args.out_edge_name, sep="\t", index=None)
    node_df.to_csv(args.out_node_name, sep="\t")

    # TODO(lev): add the size of the df to the output name


def melt_gct(gct):
    """ Melt the data_df of a connectivity gct.

    Args
    ----
    gct (GCToo object)

    Returns
    -------
    melted_df (pandas df): size = n x 3

    """
    # TODO: what about cell_id??

    assert gct.data_df.shape[0] == gct.data_df.shape[1], (
        "gct.data_df must be square. gct.data_df.shape: {}").format(gct.data_df.shape)

    # Rename the index name to 'target' before resetting it
    gct.data_df.index.name = "target"
    data_df = gct.data_df.reset_index()

    # Melt the df
    melted_df = pd.melt(data_df, id_vars="target", var_name="query")

    # Reorder the columns
    melted_df = melted_df[["query", "target", "value"]]

    return melted_df


def remove_self_cxns(melted_df_w_self_conns):

    # Remove self-connectivities
    bool_array_to_excl_self_conn = melted_df_w_self_conns.loc[:, "query"] != melted_df_w_self_conns.loc[:, "target"]
    melted_df_wo_self_conns = melted_df_w_self_conns.loc[bool_array_to_excl_self_conn, :]

    return melted_df_wo_self_conns


def network_query_of_df(melted_df, col_metadata_df, my_query, threshold, annot_fields):
    """ Main method. Works on a melted dataframe, returns node and edge dfs
    that can written to files and used by Cytoscape and other network
    visualizations. Column metadata from the GCT is used for annotating
    the nodes with annot_fields.

    For each query, all first-order connections above the threshold are found.
    Then, all connections between those elements are found; these are called
    the second-order connections.

    N is the total number of connections, or edges -- both first and
    second-order. n is the number of nodes.

    Args
    ----
    melted_df (pandas df)
    col_metadata_df (pandas df)
    my_query (list of strings): pert_ids
    threshold (float)
    annot_fields (list of strings): metadata headers

    Returns
    -------
    node_df (pandas df): size = n x a+1, where a is the length of annot_fields
    edge_df (pandas df): size = N x 5 ["query", "target", "value", "abs_value", "sign"]

    """

    # If no query provided, get all connections above threshold
    if my_query is None:
        edge_df, node_df = get_all_cxns_above_thresh(melted_df, col_metadata_df, threshold, annot_fields)

    elif type(my_query) is list:

        # Initialize lists
        list_of_edge_dfs = []
        list_of_queries = []
        all_first_order_targets = []

        # For each query, get its first-order connections
        for query_pert_iname in my_query:

            # Convert pert_iname to pert_id, append to list
            query_pert_id = convert_pert_iname_to_pert_id(query_pert_iname, col_metadata_df)
            list_of_queries.append(query_pert_id)

            # Get first-order connections
            first_order_cxns_df, first_order_targets = get_cxns_to_query_above_thresh(melted_df, query_pert_id, threshold)

            # Append targets and dfs to lists
            all_first_order_targets += first_order_targets
            list_of_edge_dfs.append(first_order_cxns_df)

        # Get second-order connections between all first-order connections
        second_order_cxns_df = get_second_order_cxns(melted_df, all_first_order_targets, threshold)
        list_of_edge_dfs.append(second_order_cxns_df)

        # Concatenate first and second-order connections & remove duplicate rows
        edge_df = create_edge_df(list_of_edge_dfs)

        # Create node_df
        all_nodes_not_unique = list_of_queries + all_first_order_targets
        all_nodes = list(set(all_nodes_not_unique))
        node_df = create_node_df(col_metadata_df, all_nodes, annot_fields)

    else:
        msg = "my_query must be a list. my_query: {}".format(my_query)
        logger.error(msg)
        raise Exception(msg)

    return node_df, edge_df


def get_all_cxns_above_thresh(melted_df, col_meta_df, threshold, annot_fields):

    # Create edge_df
    edge_df = melted_df.loc[abs(melted_df.loc[:, "value"]) >= threshold, :]

    # Add two columns: abs_value and sign
    edge_df.loc[:, "abs_value"] = abs(edge_df.loc[:, "value"]).values
    edge_df.loc[:, "sign"] = np.sign(edge_df.loc[:, "value"].values)

    # Extract nodes
    nodes_not_unique = edge_df.loc[:, "query"].tolist() + edge_df.loc[:, "target"].tolist()
    nodes = list(set(nodes_not_unique))

    # Create node_df
    node_df = create_node_df(col_meta_df, nodes, annot_fields)

    return edge_df, node_df


def convert_pert_iname_to_pert_id(query_pert_iname, col_metadata_df):

    query_id = col_metadata_df.index[col_metadata_df.loc[:, "pert_iname"] == query_pert_iname]

    assert query_id.size > 0, (
        "query_pert_iname {} not found in column metadata.").format(query_pert_iname)

    assert query_id.size == 1, (
        "More than 1 pert_id corresponds to pert_iname {}. " +
        "Mapping should be 1-to-1. query_id.values: {}").format(
        query_pert_iname, query_id.values)

    return query_id[0]


def get_cxns_to_query_above_thresh(melted_df, my_query, threshold):
    """ Get first order connections for a single query.

    Args
    ----
    melted_df (pandas df)
    my_query (string)
    threshold (float)

    Returns
    -------
    first_order_cxns_df (pandas df)
    first_order_targets (list)

    """
    # Make sure query is valid
    assert my_query in melted_df.loc[:, "query"].values, (
        "Query {} not present".format(my_query))

    # Get first order connections
    bool_array_for_my_query = (melted_df.loc[:, "query"] == my_query) & (abs(melted_df.loc[:, "value"]) >= threshold)
    first_order_cxns_df = melted_df.loc[bool_array_for_my_query, :]

    # Return first order targets as a list
    first_order_targets = first_order_cxns_df.loc[:, "target"].tolist()

    return first_order_cxns_df, first_order_targets


def create_node_df(col_meta_df, nodes, annot_fields):

    node_df = col_meta_df.loc[nodes, annot_fields]
    node_df.index.name = "pert_id"

    return node_df


def get_second_order_cxns(melted_df, first_order_targets, threshold):
    """ Get second-order connections above the threshold, i.e. only connections
    between first_order_targets.

    Args
    ----
    melted_df (pandas df)
    first_order_targets (list)
    threshold (float)

    Returns
    -------
    second_order_cxns_df (pandas df)

    """
    # Initialize a list of all False
    bool_array_of_strong_secondary_cxns = pd.Series(False, index=melted_df.index)

    # Loop over all possible permutations between these targets
    for target1, target2 in itertools.permutations(first_order_targets, 2):

        # Get boolean array of the second order connections
        this_bool_array_of_secondary_cxns = (melted_df.loc[:, "query"] == target1) & (melted_df.loc[:, "target"] == target2)

        # Only return those that are above the threshold
        this_bool_array_of_strong_secondary_cxns = this_bool_array_of_secondary_cxns & (abs(melted_df.loc[:, "value"]) >= threshold)

        # Use OR operation on this bool_array and the one external to the for-loop
        bool_array_of_strong_secondary_cxns = bool_array_of_strong_secondary_cxns | this_bool_array_of_strong_secondary_cxns

    # Return second order connections
    second_order_cxns_df = melted_df.loc[bool_array_of_strong_secondary_cxns, :]
    return second_order_cxns_df


def create_edge_df(list_of_edge_dfs):
    """ Concatenate first and second-order connections to form edge_df.
    Remove duplicate rows, and add columns for absolute_value and sign.

    Args
    ----
    list_of_edge_dfs (list of pandas df)

    Returns
    -------
    edge_df (pandas df): each row is a connection

    """
    # Concatenate
    edge_df_w_dup = pd.concat(list_of_edge_dfs)

    # Remove duplicates
    edge_df = edge_df_w_dup.drop_duplicates()
    if edge_df_w_dup.equals(edge_df):
        logger.info("Duplicate rows removed from edge_df_w_dup.")

    # Add two columns: abs_value and sign
    edge_df.loc[:, "abs_value"] = abs(edge_df.loc[:, "value"]).values
    edge_df.loc[:, "sign"] = np.sign(edge_df.loc[:, "value"].values)

    return edge_df


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)