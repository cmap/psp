"""
sip.py

Computes connectivity (KS-test or percentile scores) between a test similarity
gct and a background similarity gct. The default output is signed connectivity,
which means that the connectivity score is artifically made negative if the
median of the test distribution is less than the median of the background
distribution.

Required inputs are paths to the test and background gct files. Output is a
connectivity gct.

Metadata for the connectivity gct comes from the test gct. The dimensions will
also be the same, except the connectivity gct will not include rows
that are not also in the the background gct. Therefore, it is important that the
rows of the background gct include the rows (i.e. targets) of the test gct;
any target that is not in the background gct will not have a background
distribution, and therefore connectivity cannot be computed for that target.

N.B. The connectivity gct results will be sorted (case-insensitively, which is
the Python default).

"""

import logging
import argparse
import sys
import pandas as pd
import numpy as np
from scipy import stats

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.write_gct as wg

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

# Set up logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

CONNECTIVITY_METRIC_FIELD = "connectivity_metric"
QUERY_FIELD_NAME = "query_field"
TARGET_FIELD_NAME = "target_field"


def build_parser():
    """ Build argument parser. """

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--test_gct_path", "-t", required=True,
                        help="path to input gct file")
    parser.add_argument("--bg_gct_path", "-b", required=True,
                        help="path to background gct file")

    # Optional args
    parser.add_argument("--out_name", "-o", default="sip_output.gct",
                        help="what to name the output connectivity file")
    parser.add_argument("--connectivity_metric", "-c", default="ks_test",
                        choices=["ks_test", "percentile_score"],
                        help="metric to use for computing connectivity")
    parser.add_argument("--fields_to_aggregate_in_test_gct_queries", "-tfq",
                        nargs="*", default=["pert_id", "cell_id", "pert_time"],
                        help="list of metadata fields in the columns of the test gct to aggregate")
    parser.add_argument("--fields_to_aggregate_in_test_gct_targets", "-tft",
                        nargs="*", default=["pert_id", "cell_id", "pert_time"],
                        help="list of metadata fields in the rows of the test gct to aggregate")
    parser.add_argument("--fields_to_aggregate_in_bg_gct", "-bf",
                        nargs="*", default=["pert_id", "cell_id", "pert_time"],
                        help="list of metadata fields in the bg gct to aggregate")
    parser.add_argument("--separator", "-s", type=str, default=":",
                        help="string separator for aggregating fields together")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):
    """ The main method. """

    # Read test gct
    test_gct = parse(args.test_gct_path, convert_neg_666=False, make_multiindex=True)

    # Read bg_gct
    bg_gct = parse(args.bg_gct_path, convert_neg_666=False, make_multiindex=True)

    # Create an aggregated metadata field for index and columns of both gcts
    # and sort by that field
    (test_df, bg_df) = prepare_multi_index_dfs(
        test_gct.multi_index_df, bg_gct.multi_index_df,
        args.fields_to_aggregate_in_test_gct_queries,
        args.fields_to_aggregate_in_test_gct_targets,
        args.fields_to_aggregate_in_bg_gct,
        QUERY_FIELD_NAME,
        TARGET_FIELD_NAME,
        args.separator)

    # Check symmetry
    (is_test_df_sym, _) = check_symmetry(test_gct.multi_index_df, bg_gct.multi_index_df)

    # Compute connectivity
    (conn_mi_df, signed_conn_mi_df) = compute_connectivities(
        test_df, bg_df, QUERY_FIELD_NAME, TARGET_FIELD_NAME, TARGET_FIELD_NAME,
        args.connectivity_metric, is_test_df_sym)

    # Convert multi-index to component dfs in order to write output gct
    (signed_data_df, signed_row_metadata_df, signed_col_metadata_df) = (
        GCToo.multi_index_df_to_component_dfs(
            signed_conn_mi_df, rid=TARGET_FIELD_NAME, cid=QUERY_FIELD_NAME))

    # Append to queries a new column saying what connectivity metric was used
    add_connectivity_metric_to_metadata(signed_col_metadata_df, args.connectivity_metric, CONNECTIVITY_METRIC_FIELD)
    add_connectivity_metric_to_metadata(signed_row_metadata_df, args.connectivity_metric, CONNECTIVITY_METRIC_FIELD)

    # Create gct and write it to file
    conn_gct = GCToo.GCToo(data_df=signed_data_df, row_metadata_df=signed_row_metadata_df, col_metadata_df=signed_col_metadata_df)
    wg.write(conn_gct, args.out_name, data_null="NaN", filler_null="NaN", metadata_null="NaN")


def prepare_multi_index_dfs(test_df, bg_df, fields_to_aggregate_in_test_gct_queries,
                            fields_to_aggregate_in_test_gct_targets,
                            fields_to_aggregate_in_bg_gct, query_field_name,
                            target_field_name, sep):
    """
    This functions prepare the test and background multi-index dfs before
    connectivity can be computed.

    For each fields_to_aggregate_* argument, this function adds a level to the
    corresponding multi_index that is an aggregation of the entries in
    fields_to_aggregate_*. The name of this level (aggregated_field_name) is
    returned as out_*_field. If fields_to_aggregate_* has only 1 element,
    a new field is still made but it's just a copy of what's in
    fields_to_aggregate_*[0].

    This fcn also sorts the dfs by out_field (sorting is necessary for efficient
    indexing later).

    Args:
        test_df (multi-index pandas df)
        bg_df (multi-index pandas df)
        fields_to_aggregate_in_test_gct_queries (list of strings)
        fields_to_aggregate_in_test_gct_targets (list of strings)
        fields_to_aggregate_in_bg_gct (list of strings)
        query_field_name (string)
        target_field_name (string)
        sep (string)

    Returns:
        out_test_df (multi-index pandas df)
        out_bg_df (multi-index pandas df)

    """
    # Check if columns are multi-index or regular index
    if isinstance(test_df.columns, pd.core.index.MultiIndex):

        # If no aggregation fields were provided, use cid
        if len(fields_to_aggregate_in_test_gct_queries) == 0:
            fields_to_aggregate_in_test_gct_queries = ["cid"]

        # Create aggregated level in test_df for queries (columns)
        (_, test_df_columns) = add_aggregated_level_to_multi_index(
            test_df.columns, fields_to_aggregate_in_test_gct_queries, query_field_name, sep)

    else:
        logger.info(("No metadata provided for test GCT columns. " +
                     "Using ids as perturbation identifiers."))
        test_df_columns = index_to_multi_index(test_df.columns, query_field_name)

    # Check if index is multi-index or regular index
    if isinstance(test_df.index, pd.core.index.MultiIndex):

        # If no aggregation fields were provided, use rid
        if len(fields_to_aggregate_in_test_gct_targets) == 0:
            fields_to_aggregate_in_test_gct_targets = ["rid"]

        # Create aggregated level in test_df for targets (rows)
        (_, test_df_index) = add_aggregated_level_to_multi_index(
            test_df.index, fields_to_aggregate_in_test_gct_targets, target_field_name, sep)

    else:
        logger.info(("No metadata provided for test GCT index. " +
                     "Using ids as perturbation identifiers."))
        test_df_index = index_to_multi_index(test_df.index, target_field_name)

    # Create out_test_df
    out_test_df = pd.DataFrame(test_df.values, index=test_df_index, columns=test_df_columns)

    # Sort out_test_df
    out_test_df.sort_index(level=target_field_name, axis=0, inplace=True)
    out_test_df.sort_index(level=query_field_name, axis=1, inplace=True)

    # Check if columns are multi-index or regular index
    if isinstance(bg_df.columns, pd.core.index.MultiIndex):

        # If no aggregation fields were provided, use cid
        if len(fields_to_aggregate_in_bg_gct) == 0:
            (_, bg_df_columns) = add_aggregated_level_to_multi_index(
                bg_df.columns, ["cid"], target_field_name, sep)

        else:
            # Create aggregated level in bg_df
            (_, bg_df_columns) = add_aggregated_level_to_multi_index(
                bg_df.columns, fields_to_aggregate_in_bg_gct, target_field_name, sep)

    else:
        logger.info(("No metadata provided for background GCT columns. " +
                     "Using ids as perturbation identifiers."))
        bg_df_columns = index_to_multi_index(bg_df.columns, target_field_name)

    # Check if index multi-index or regular index
    if isinstance(bg_df.index, pd.core.index.MultiIndex):

        # If no aggregation fields were provided, use rid
        if len(fields_to_aggregate_in_bg_gct) == 0:
            (_, bg_df_index) = add_aggregated_level_to_multi_index(
                bg_df.index, ["rid"], target_field_name, sep)

        else:
            # Create aggregated level in bg_df
            (_, bg_df_index) = add_aggregated_level_to_multi_index(
                bg_df.index, fields_to_aggregate_in_bg_gct, target_field_name, sep)

    else:
        logger.info(("No metadata provided for background GCT index. " +
                     "Using ids as perturbation identifiers."))
        bg_df_index = index_to_multi_index(bg_df.index, target_field_name)

    # Create out_bg_df
    out_bg_df = pd.DataFrame(bg_df.values, index=bg_df_index, columns=bg_df_columns)

    # Sort out_bg_df
    out_bg_df.sort_index(level=target_field_name, axis=0, inplace=True)
    out_bg_df.sort_index(level=target_field_name, axis=1, inplace=True)

    return out_test_df, out_bg_df


def check_symmetry(test_mi_df, bg_mi_df):
    """ Checks if test multi-index df and bg multi-index dfs are symmetric.
    If so, extraction of values later has to be different to avoid double-counting.

    Args:
        test_mi_df (pandas multi-index df)
        bg_mi_df (pandas multi-index df)

    Returns:
        is_test_df_sym (bool)
        is_bg_df_sym (bool)

    """
    # Check if test_mi_df is symmetric
    if test_mi_df.index.equals(test_mi_df.columns):
        is_test_df_sym = True
        logger.info("test_mi_df has the same index and columns, so it will be considered symmetric.")
    else:
        is_test_df_sym = False

    # Check if bg_mi_df is symmetric
    if bg_mi_df.index.equals(bg_mi_df.columns):
        is_bg_df_sym = True
        logger.info("bg_mi_df has the same index and columns, so it will be considered symmetric.")
    else:
        is_bg_df_sym = False

    # TODO(lev): should be able to support a non-symmetric background matrix
    # TODO(lev): when we do this, we should also separate bg_gct_field into bg_gct_query_field and bg_gct_target_field
    assert is_bg_df_sym, "bg_mi_df must be symmetric!"

    return is_test_df_sym, is_bg_df_sym


def compute_connectivities(test_df, bg_df, test_gct_query_field, test_gct_target_field, bg_gct_field, connectivity_metric, is_test_df_sym):
    """ Compute all connectivities for a single test_df and a single bg_df.

    Args:
        test_df (multi-index or regular pandas df): m x n, where n is the # of queries, m is the # of targets
        bg_df (multi-index or regular pandas df): M x M, where M includes m
        test_gct_query_field (string)
        test_gct_target_field (string)
        bg_gct_field (string)
        connectivity_metric (string)
        is_test_df_sym (bool)

    Returns:
        conn_df (multi-index or regular pandas df): m x n, where n is the # of queries, m is the # of targets
        signed_conn_df (multi-index or regular  pandas df): m x n, where n is the # of queries, m is the # of targets

    """
    logger.info("Computing connectivities...")

    # Extract queries from test_df columns and targets from test_df index
    queries = test_df.columns.get_level_values(test_gct_query_field).unique()

    # Extract targets from both test_df and bg_df
    bg_targets = bg_df.index.get_level_values(bg_gct_field).unique()
    test_targets = test_df.index.get_level_values(test_gct_target_field).unique()

    # Get intersection of targets
    targets = np.intersect1d(test_targets, bg_targets)
    assert targets.size > 0, ("There are no targets in common between the test and bg dfs.\n" +
        "test_targets: {}, bg_targets: {}".format(test_targets, bg_targets))

    logger.info("{} queries, {} targets".format(len(queries), len(targets)))

    # Column multi-index of output df is the unique version of test_df.columns
    conn_df_columns = pd.MultiIndex.from_tuples(np.unique(test_df.columns.values))
    conn_df_columns.names = test_df.columns.names

    # Sort by test_gct_query_field to ensure that it still aligns with queries
    sorted_conn_df_columns = conn_df_columns.sortlevel(test_gct_query_field)[0]
    assert np.array_equal(
        sorted_conn_df_columns.get_level_values(test_gct_query_field), queries), (
        "queries do not match the {test_gct_query_field} level in sorted_conn_df_columns. " +
        "queries:\n{queries}\nsorted_conn_df_columns.get_level_values({test_gct_query_field}):\n{sorted_conn_df}").format(
        test_gct_query_field=test_gct_query_field,
        queries=queries,
        sorted_conn_df=sorted_conn_df_columns.get_level_values(test_gct_query_field))

    # Row multi-index of output df is the unique version of test_df.index, but
    # only for targets also in bg_df
    conn_df_index = pd.MultiIndex.from_tuples(np.unique(test_df.index[test_df.index.get_level_values(test_gct_target_field).isin(targets)].values))
    conn_df_index.names = test_df.index.names

    # Sort by test_gct_target_field to ensure that it still aligns with targets
    sorted_conn_df_index = conn_df_index.sortlevel(test_gct_target_field)[0]
    assert np.array_equal(
        sorted_conn_df_index.get_level_values(test_gct_target_field), targets), (
        "targets do not match the {test_gct_target_field} level in sorted_conn_df_index. " +
        "targets:\n{targets}\nsorted_conn_df_index.get_level_values({test_gct_target_field}):\n{sorted_conn_df}").format(
        test_gct_target_field=test_gct_target_field,
        targets=targets,
        sorted_conn_df=sorted_conn_df_index.get_level_values(test_gct_target_field))

    # Initialize conn_df and signed_conn_df :: len(targets) x len(queries)
    # Use regular rather than multi-index dfs to make insertion of connectivity
    # values easier; convert to multi-indices after the for-loop
    conn_df = pd.DataFrame(np.zeros([len(targets), len(queries)]) * np.nan,
                           index=targets, columns=queries)
    signed_conn_df = conn_df.copy()

    for target in targets:

        # Extract background values
        bg_vals = extract_bg_vals_from_sym(target, bg_gct_field, bg_df)

        # Make sure bg_vals has at least 1 element before continuing
        if len(bg_vals) > 0:

            for query in queries:
                logger.debug("query: {}, target: {}".format(query, target))

                # Extract test values
                test_vals = extract_test_vals(query, target, test_gct_query_field,
                                              test_gct_target_field, test_df, is_test_df_sym)

                # Make sure test_vals has at least 1 element before continuing
                if len(test_vals) > 0:

                    if connectivity_metric == "ks_test":

                        # Compute single connectivity
                        (ks_stat, _) = ks_test_single(test_vals, bg_vals)
                        conn_df.loc[target, query] = ks_stat

                        # Compute signed connectivity as well
                        signed_ks_stat = add_sign_to_conn(ks_stat, test_vals, bg_vals)
                        signed_conn_df.loc[target, query] = signed_ks_stat

                    elif connectivity_metric == "percentile_score":

                        # Compute single connectivity
                        conn = percentile_score_single(test_vals, bg_vals)
                        conn_df.loc[target, query] = conn

                        # Compute signed connectivity as well
                        signed_conn = add_sign_to_conn(conn, test_vals, bg_vals)
                        signed_conn_df.loc[target, query] = signed_conn

                    else:
                        err_msg = ("connectivity metric must be either ks_test or " +
                                   "percentile_score. connectivity_metric: {}").format(connectivity_metric)
                        logger.error(err_msg)
                        raise(Exception(err_msg))

    # Replace the regular indices with multi-indices
    conn_df.index = sorted_conn_df_index
    conn_df.columns = sorted_conn_df_columns
    signed_conn_df.index = sorted_conn_df_index
    signed_conn_df.columns = sorted_conn_df_columns

    return conn_df, signed_conn_df


def ks_test_single(test_vals, bg_vals):
    """ Compute KS-test statistic for one pair of test values and background
    values.

    Args:
        test_vals (numpy array)
        bg_vals (numpy array)

    Returns:
        ks_stat (float)
        pval (float)

    """
    # Do KS-test
    try:
        (ks_stat, pval) = stats.ks_2samp(test_vals, bg_vals)

    # Return NaN if test fails
    except ValueError:
        logger.warning("KS-test failed.")
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


def extract_test_vals(query, target, query_field, target_field, test_df, is_test_df_sym):
    """ Extract values that has query in the columns and target in the rows.

    Args:
        query (string)
        target (string)
        query_field (string): name of multiindex level in which to find query
        target_field (string): name of multiindex level in which to find target
        test_df (pandas multi-index df)
        is_test_df_sym (bool): only matters if query == target; set to True to
            avoid double-counting in the case of a symmetric matrix

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

    # If query == target AND the matrix is symmetric, need to take only triu
    # of the extracted values in order to avoid double-counting
    if query == target and is_test_df_sym:
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


def index_to_multi_index(regular_index, new_level_name):
    out_mi = pd.MultiIndex.from_arrays([regular_index.values, regular_index.values],
                                       names=["id", new_level_name])

    return out_mi


def add_aggregated_level_to_multi_index(mi, levels_to_aggregate, aggregated_level_name, sep):
    """ Create a new level by aggregating values in other multi-index levels.

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
    # N.B. Convert each level to a string in order to produce a string ultimately
    list_of_levels = [mi.get_level_values(level).values.astype(str) for level in levels_to_aggregate]

    # Zip each of the levels together into a tuple
    list_of_aggregated_tuples = zip(*list_of_levels)

    # Join the tuple strings to make the aggregated strings
    aggregated_strings = [sep.join(agg_tuple) for agg_tuple in list_of_aggregated_tuples]

    # Convert each level of multi-index to its own tuple
    levels_as_tuples = zip(*mi.ravel())

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


def add_connectivity_metric_to_metadata(metadata_df, connectivity_metric, connectivity_metric_field):
    """ Append a connectivity_metric_field to metadata_df that is the
    connectivity metric repeated.

    Args:
        metadata_df (pandas df)
        connectivity_metric (string)
        connectivity_metric_field (string): name of the connectivity_metric column

    Returns:
        metadata_df (pandas df): modified in-place

    """
    metadata_df.loc[:, connectivity_metric_field] = connectivity_metric

    return None


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)