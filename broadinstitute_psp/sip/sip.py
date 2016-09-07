"""
sip.py

Computes connectivity (KS-test or percentile scores) between a test similarity
gct and a background similarity gct. The default output is signed connectivity,
which means that the connectivity score is artifically made negative if the
median of the test distribution is less than the median of the background
distribution.

Required inputs are paths to the test and background gct files. Output is a
connectivity gct.

The dimensions of the connectivity gct will be equal to the dimensions of the
test gct. It is important that the rows of the background gct include the rows
(i.e. targets) of the test gct; any target that is not in the background gct
will not have a background distribution, and therefore connectivity cannot be
computed for that target.

N.B. The connectivity gct results will be sorted (case-insensitively, which is
the Python default).

"""

import logging
import argparse
import sys
import pandas as pd
import numpy as np
from scipy import stats

import utils.setup_logger as setup_logger
import utils.psp_utils as utils
import cmap.io.GCToo.GCToo as GCToo
import cmap.io.GCToo.write_gctoo as wg

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

# Set up logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

CONNECTIVITY_METRIC_FIELD = "connectivity_metric"

def build_parser():
    """ Build argument parser. """

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--test_gct_path", "-t", required=True,
                        help="path to input gct file")
    parser.add_argument("--bg_gct_path", "-b", required=True,
                        help="path to output directory")

    # Optional args
    parser.add_argument("--out_name", "-o", default="sip_output.gct",
                        help="what to name the output connectivity file")
    parser.add_argument("--connectivity_metric", "-c", default="ks_test",
                        choices=["ks_test", "percentile_score"],
                        help="metric to use for computing connectivity")
    parser.add_argument("--psp_config_path", "-p",
                        default="~/psp_production.cfg",
                        help="filepath to PSP config file")
    parser.add_argument("--fields_to_aggregate_in_test_gct_queries", "-tfq",
                        nargs="+", default=["pert_iname"],
                        help="list of metadata fields in the columns of the test gct to aggregate")
    parser.add_argument("--fields_to_aggregate_in_test_gct_targets", "-tft",
                        nargs="+", default=["pert_iname"],
                        help="list of metadata fields in the rows of the test gct to aggregate")
    parser.add_argument("--fields_to_aggregate_in_bg_gct", "-bf",
                        nargs="+", default=["pert_iname"],
                        help="list of metadata fields in the bg gct to aggregate")
    parser.add_argument("--aggregated_field_name", "-a",
                        type=str, default="aggregated",
                        help="what to call the aggregated metadata field")
    parser.add_argument("--separator", "-s", type=str, default=":",
                        help="string separator for the aggregated field")
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):
    """ The main method. """

    # Read test gct and config file
    test_gct = utils.read_gct_and_config_file(args.test_gct_path, args.psp_config_path)[0]

    # Read bg_gct
    bg_gct = utils.read_gct_and_config_file(args.bg_gct_path, args.psp_config_path)[0]

    # Create an aggregated metadata field for index and columns of both gcts
    # and sort by that field
    (test_df, bg_df, test_query_field, test_target_field, bg_field) = prepare_multi_index_dfs(
        test_gct.multi_index_df, bg_gct.multi_index_df,
        args.fields_to_aggregate_in_test_gct_queries,
        args.fields_to_aggregate_in_test_gct_targets,
        args.fields_to_aggregate_in_bg_gct,
        args.aggregated_field_name,
        args.separator)

    # Compute connectivity
    (conn_df, signed_conn_df) = compute_connectivities(
        test_df, bg_df, test_query_field, test_target_field, bg_field, args.connectivity_metric)

    # Create connectivity metadata dfs (rows are targets, columns are queries)
    row_metadata_df = create_connectivity_metadata_df(
        conn_df.index.values, args.connectivity_metric, CONNECTIVITY_METRIC_FIELD,
        "target", "query_or_target")
    col_metadata_df = create_connectivity_metadata_df(
        conn_df.columns.values, args.connectivity_metric, CONNECTIVITY_METRIC_FIELD,
        "query", "query_or_target")

    # Create output gct
    conn_gct = GCToo.GCToo(data_df=signed_conn_df, row_metadata_df=row_metadata_df, col_metadata_df=col_metadata_df)

    # Write gctoo to file
    wg.write(conn_gct, args.out_name, data_null="NaN", filler_null="NaN", metadata_null="NaN")


def prepare_multi_index_dfs(test_df, bg_df, fields_to_aggregate_in_test_gct_queries,
                            fields_to_aggregate_in_test_gct_targets,
                            fields_to_aggregate_in_bg_gct, aggregated_field_name, sep):
    """
    This functions prepare the test and background multi-index dfs before
    connectivity can be computed.

    For each fields_to_aggregate_* argument, if it has more than 1 element,
    this function adds a level to the corresponding multi_index that is an
    aggregation of the entries in fields_to_aggregate_*. The name of this level
    (aggregated_field_name) is returned as out_*_field. If fields_to_aggregate_*
    has only 1 element, this single element is returned as out_*_field.

    This fcn also sorts the dfs by out_field (sorting is necessary for efficient
    indexing later).

    Args:
        test_df (multi-index pandas df)
        bg_df (multi-index pandas df)
        fields_to_aggregate_in_test_gct_queries (list of strings)
        fields_to_aggregate_in_test_gct_targets (list of strings)
        fields_to_aggregate_in_bg_gct (list of strings)
        aggregated_field_name (string)
        sep (string)

    Returns:
        out_test_df (multi-index pandas df)
        out_bg_df (multi-index pandas df)
        out_test_query_field (string)
        out_test_target_field (string)
        out_bg_field (string)

    """
    # Create aggregated level for queries (columns) in test_df if needed
    if len(fields_to_aggregate_in_test_gct_queries) > 1:
        test_df_columns = add_aggregated_level_to_multi_index(
            test_df.columns, fields_to_aggregate_in_test_gct_queries, aggregated_field_name, sep)
        out_test_query_field = aggregated_field_name
    else:
        test_df_columns = test_df.columns
        out_test_query_field = fields_to_aggregate_in_test_gct_queries[0]

    # Sort test_df_columns
    sorted_test_df_columns = test_df_columns.sort_values()

    # Create aggregated level for targets (rows) in test_df if needed
    if len(fields_to_aggregate_in_test_gct_targets) > 1:
        test_df_index = add_aggregated_level_to_multi_index(
            test_df.index, fields_to_aggregate_in_test_gct_targets, aggregated_field_name, sep)
        out_test_target_field = aggregated_field_name
    else:
        test_df_index = test_df.index
        out_test_target_field = fields_to_aggregate_in_test_gct_targets[0]

    # Sort test_df_index
    sorted_test_df_index = test_df_index.sort_values()

    # Create out_test_df
    out_test_df = pd.DataFrame(test_df.values, index=sorted_test_df_index, columns=sorted_test_df_columns)

    # Create aggregated level in bg_df if needed
    if len(fields_to_aggregate_in_bg_gct) > 1:
        bg_df_columns = add_aggregated_level_to_multi_index(
            bg_df.columns, fields_to_aggregate_in_bg_gct, aggregated_field_name, sep)
        bg_df_index = add_aggregated_level_to_multi_index(
            bg_df.index, fields_to_aggregate_in_bg_gct, aggregated_field_name, sep)
        out_bg_field = aggregated_field_name
    else:
        bg_df_columns = bg_df.columns
        bg_df_index = bg_df.index
        out_bg_field = fields_to_aggregate_in_bg_gct[0]

    # Sort bg_df_columns and bg_df_index
    sorted_bg_df_columns = bg_df_columns.sort_values()
    sorted_bg_df_index = bg_df_index.sort_values()

    # Create out_bg_df
    out_bg_df = pd.DataFrame(bg_df.values, index=sorted_bg_df_index, columns=sorted_bg_df_columns)

    return out_test_df, out_bg_df, out_test_query_field, out_test_target_field, out_bg_field


def compute_connectivities(test_df, bg_df, test_gct_query_field, test_gct_target_field, bg_gct_field, connectivity_metric):
    """ Compute all connectivities for a single test_df and a single bg_df.

    Args:
        test_df (multi-index pandas df): m x n, where n is the # of queries, m is the # of targets
        bg_df (multi-index pandas df): M x M, where M includes m entries
        test_gct_query_field (string)
        test_gct_target_field (string)
        bg_gct_field (string)
        connectivity_metric (string)

    Returns:
        conn_df (pandas df): m x n, where n is the # of queries, m is the # of targets

    """
    logger.info("Computing connectivities...")

    # TODO(lev): should be able to support a non-symmetric background matrix
    # TODO(lev): when we do this, we should also separate bg_gct_field into bg_gct_query_field and bg_gct_target_field
    is_sym = bg_df.index.equals(bg_df.columns)
    assert is_sym, "bg_df must be symmetric"

    # Extract queries from test_df columns and targets from test_df index
    queries = test_df.columns.get_level_values(test_gct_query_field).unique()
    targets = test_df.index.get_level_values(test_gct_target_field).unique()
    logger.info("{} queries, {} targets".format(len(queries), len(targets)))

    # Initialize conn_df :: len(targets) x len(queries)
    conn_df = pd.DataFrame(np.zeros([len(targets), len(queries)]) * np.nan,
                           index=targets, columns=queries)
    signed_conn_df = conn_df.copy()

    for target in targets:
        bg_vals = extract_bg_vals_from_sym(target, bg_gct_field, bg_df)
        for query in queries:
            logger.debug("query: {}, target: {}".format(query, target))
            test_vals = extract_test_vals(query, target, test_gct_query_field, test_gct_target_field, test_df)
        
            if connectivity_metric == "ks_test":

                # Compute single connectivity
                (ks_stat, pval) = ks_test_single(test_vals, bg_vals)
                conn_df.loc[target, query] = ks_stat

                # TODO(lev): figure out what to do with pvals

                # Compute signed connectivity as well
                signed_ks_stat = add_sign_to_conn(ks_stat, test_vals, bg_vals)
                signed_conn_df.loc[target, query] = signed_ks_stat

            elif connectivity_metric == "percentile_score":

                # Compute single connectivity
                conn = percentile_score_single(test_vals, bg_vals)
                conn_df.loc[target, query] = conn

            else:
                err_msg = ("connectivity metric must be either ks_test or " +
                           "percentile_score. connectivity_metric: {}").format(connectivity_metric)
                logger.error(err_msg)
                raise(Exception(err_msg))

    return conn_df, signed_conn_df


def ks_test_single(test_vals, bg_vals, min_number_of_elements=2):
    """
    Compute KS-test statistic for one pair of test values and background values.

    min_number_of_elements can be used to make sure that each distribution has
    enough elements for the result of the KS-test to be meaningful.
    """

    # Check that each distribution has some minimum number of elements
    if len(test_vals) >= min_number_of_elements and len(bg_vals) >= min_number_of_elements:

        # Do KS-test
        try:
            (ks_stat, pval) = stats.ks_2samp(test_vals, bg_vals)

        # Return NaN if test fails for some reason
        except ValueError as e:
            logger.warning("KS-test failed.")
            ks_stat = np.nan
            pval = np.nan

    else:
        ks_stat = np.nan
        pval = np.nan

    return ks_stat, pval


def percentile_score_single(test_vals, bg_vals):
    """ For each value in test_vals, compute its percentile score compared
    to bg_vals.

    Args:
        test_vals (numpy array)
        bg_vals (numpy array)

    Returns:
        out_score (float)

    """

    # Compute percentile score for each value in test_vals
    percentile_scores = [stats.percentileofscore(bg_vals, test_val, kind="rank") for test_val in test_vals]

    # Take mean of percentile scores
    out_score = np.mean(percentile_scores)

    return out_score


def add_sign_to_conn(conn, test_vals, bg_vals):
    """
    If median of test_vals is less than the median of bg_vals,
    return the signed connectivity, i.e. signed_conn = conn * -1.
    """

    if (np.median(test_vals) - np.median(bg_vals)) >= 0:
        signed_conn = conn
    else:
        signed_conn = conn * -1

    return signed_conn


def extract_test_vals(query, target, query_field, target_field, test_df):
    """ Extract values that has query in the columns and target in the rows.

    Args:
        query (string)
        target (string)
        query_field (string): name of multiindex level in which to find query
        target_field (string): name of multiindex level in which to find target
        test_df (pandas multi-index df)

    Returns:
        vals (numpy array)

    """
    assert query in test_df.columns.get_level_values(query_field), (
        "query {} is not in the {} level of the columns of test_df.".format(
            query, query_field))

    assert target in test_df.index.get_level_values(target_field), (
        "target {} is not in the {} level of the index of test_df.".format(
            target, target_field))

    # Extract elements where query is in columns and target is in rows
    target_in_rows_query_in_cols_df = test_df.loc[
            test_df.index.get_level_values(target_field) == target,
            test_df.columns.get_level_values(query_field) == query]

    # If query == target, need to take only triu of the extracted values in
    # order to avoid double-counting
    if query == target:
        mask = np.triu(np.ones(target_in_rows_query_in_cols_df.shape), k=1).astype(np.bool)
        vals_with_nans = target_in_rows_query_in_cols_df.where(mask).values.flatten()
        vals = vals_with_nans[~np.isnan(vals_with_nans)]

    else:
        vals = target_in_rows_query_in_cols_df.values.flatten()

    return vals


def extract_bg_vals_from_sym(target, multi_index_level_name, bg_df):
    """ Extract all values that have some interaction with target.

    Diagonal and lower-right triangle are excluded.

    Args:
        target (string)
        multi_index_level_name (string)
        bg_df (multi-index pandas df)

    Returns:
        vals (numpy array)

    """
    assert target in bg_df.index.get_level_values(multi_index_level_name), (
        "target {} is not in the {} level of the index of bg_df.".format(
            target, multi_index_level_name))

    assert target in bg_df.columns.get_level_values(multi_index_level_name), (
        "target {} is not in the {} level of the columns of bg_df.".format(
            target, multi_index_level_name))

    # Get values for multi_index_level_name in the index
    index_names = bg_df.index.get_level_values(multi_index_level_name)

    # Get values for multi_index_level_name in the columns
    column_names = bg_df.columns.get_level_values(multi_index_level_name)

    # Initialize list of idxs for the values we want to extract
    idxs = []

    # Loop over the rows
    for row_idx in range(len(bg_df)):
        for col_idx in range(len(bg_df)):
            # Extract element if it is in the upper-right triangle, excluding
            # the diagonal, and it contains an interaction with target
            if row_idx < col_idx and (index_names[row_idx] == target or column_names[col_idx] == target):
                idxs.append([row_idx, col_idx])

    vals = bg_df.values[tuple(np.transpose(idxs))]

    return vals

def extract_bg_vals_from_non_sym(target, multi_index_level_name, bg_df):
    """ Extract all values that have some interaction with target.

    Just get all values where target is in the rows.

    Args:
        target (string)
        multi_index_level_name (string)
        bg_df (multi-index pandas df)

    Returns:
        vals (numpy array)

    """
    assert target in bg_df.index.get_level_values(multi_index_level_name), (
        "target {} is not in the {} level of the index of bg_df.".format(
            target, multi_index_level_name))

    vals = bg_df.loc[bg_df.index.get_level_values(multi_index_level_name) == target, :].values.flatten()

    return vals


def add_aggregated_level_to_multi_index(mi, levels_to_aggregate, aggregated_level_name, sep):
    """Create a new level by aggregating values in other multi-index levels.

    In addition to returning the original multi-index with the aggregated level
    added, this function also returns subset_mi, which is the subset of the
    full output that only has levels_to_aggregate and aggregated_level_name.

    Args:
        mi (multi-index for a pandas df)
        levels_to_aggregate (list or numpy array of strings)
        aggregated_level_name (string)
        sep (string): separator to use in creating aggregated strings

    Returns:
        updated_mi (multi-index for a pandas df)
        subset_mi (multi-index for a pandas df)

    """
    # Check that each level in levels_to_aggregate is in mi
    for level in levels_to_aggregate:
        assert level in mi.names, (
                "{} is not present in the names of the multi-index.".format(level))

    # Extract each level in levels_to_aggregate
    # N.B. Convert each level to a string in order to get a string at the end
    list_of_levels = [mi.get_level_values(level).values.astype(str) for level in levels_to_aggregate]

    # Zip each of the levels together into a tuple
    list_of_aggregated_tuples = zip(*list_of_levels)

    # Join the tuple strings to make the aggregated strings
    aggregated_strings = [sep.join(agg_tuple) for agg_tuple in list_of_aggregated_tuples]

    # Convert each level of multi-index to its own tuple
    levels_as_tuples= zip(*mi.ravel())

    # Add new level as a tuple
    levels_as_tuples.append(tuple(aggregated_strings))

    # Make a new multi-index with the new level
    updated_mi = pd.MultiIndex.from_tuples(
        zip(*levels_as_tuples), names=mi.names + [aggregated_level_name])

    # Create subset_mi by dropping undesired levels
    levels_to_keep = levels_to_aggregate + [aggregated_level_name]
    levels_to_drop = [level_name for level_name in updated_mi.names if level_name not in levels_to_keep]
    subset_mi = updated_mi.copy()
    for level_to_drop in levels_to_drop:
        subset_mi = subset_mi.droplevel(level_to_drop)

    return updated_mi, subset_mi


def create_connectivity_metadata_df(perts, connectivity_metric, connectivity_metric_field, query_or_target, query_or_target_field):
    """ Create a df where the first column is which connectivity metric was
    used (same value for all rows) and the second columns indicates if this
    row was a query or target (same value for all rows).

    Args:
        perts (numpy array of strings)
        connectivity_metric (string)
        connectivity_metric_field (string): name of the connectivity_metric column
        query_or_target (string)
        query_or_target_field (string): name of the query_or_target column

    Returns:
        conn_metadata_df (pandas df)

    """
    conn_metadata_df = pd.DataFrame(
        np.tile([connectivity_metric, query_or_target], (len(perts), 1)),
        index=perts, columns=[connectivity_metric_field, query_or_target_field])

    return conn_metadata_df


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)