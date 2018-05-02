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

Any metadata present in the test gct will be retained. Since aggregation of
replicates occurs, unique metadata entries will be collapsed.

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
                        help="metadata fields in the test gct columns identifying replicates")
    parser.add_argument("--fields_to_aggregate_in_test_gct_targets", "-tft",
                        nargs="*", default=["pert_id", "cell_id", "pert_time"],
                        help="metadata fields in the test gct rows identifying replicates")
    parser.add_argument("--fields_to_aggregate_in_bg_gct", "-bf",
                        nargs="*", default=["pert_id", "cell_id", "pert_time"],
                        help="metadata fields in the background gct rows AND columns identifying replicates")
    parser.add_argument("--separator", "-s", type=str, default="|",
                        help="string separator for aggregating fields together")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):
    """ The main method. """

    # Read test gct
    test_gct = parse.parse(args.test_gct_path)

    # Read bg_gct
    bg_gct = parse.parse(args.bg_gct_path)

    # Check symmetry
    (is_test_df_sym, _) = check_symmetry(test_gct.data_df, bg_gct.data_df)

    # Create an aggregated metadata field in test and background GCTs
    # that will be used to aggregate replicates
    (test_gct, bg_gct) = create_aggregated_fields_in_GCTs(
        test_gct, bg_gct,
        args.fields_to_aggregate_in_test_gct_queries,
        args.fields_to_aggregate_in_test_gct_targets,
        args.fields_to_aggregate_in_bg_gct,
        QUERY_FIELD_NAME, TARGET_FIELD_NAME, args.separator)

    # Compute connectivity
    (conn_gct, signed_conn_gct) = compute_connectivities(
        test_gct, bg_gct, QUERY_FIELD_NAME, TARGET_FIELD_NAME, TARGET_FIELD_NAME,
        args.connectivity_metric, is_test_df_sym, args.separator)

    # Append to queries a new column saying what connectivity metric was used
    add_connectivity_metric_to_metadata(signed_conn_gct.col_metadata_df, args.connectivity_metric, CONNECTIVITY_METRIC_FIELD)
    add_connectivity_metric_to_metadata(signed_conn_gct.row_metadata_df, args.connectivity_metric, CONNECTIVITY_METRIC_FIELD)

    # Write signed result to file
    wg.write(signed_conn_gct, args.out_name, data_null="NaN", filler_null="NaN", metadata_null="NaN")


def create_aggregated_fields_in_GCTs(test_gct, bg_gct, fields_to_aggregate_in_test_gct_queries,
                                     fields_to_aggregate_in_test_gct_targets,
                                     fields_to_aggregate_in_bg_gct, query_field_name,
                                     target_field_name, sep):
    """
    This function creates new metadata fields in the GCTs indicating how
    replicates should be collapsed.

    - fields_to_aggregate_in_test_gct_queries must be in the columns of the test_gct
    - fields_to_aggregate_in_test_gct_targets must be in the rows of the test_gct
    - fields_to_aggregate_in_bg_gct must be in the rows AND columns of the bg_gct
        (currently only support symmetric background GCTs)

    If any of these arguments is an empty list, then the ids are used.

    Args:
        test_gct (GCToo)
        bg_gct (GCToo)
        fields_to_aggregate_in_test_gct_queries (list of strings)
        fields_to_aggregate_in_test_gct_targets (list of strings)
        fields_to_aggregate_in_bg_gct (list of strings)
        query_field_name (string)
        target_field_name (string)
        sep (string)

    Returns:
        test_gct (GCToo): with one row and one column metadata field appended
        bg_gct (GCToo): with one row and one column metadata field appended

    """
    # Check if we have any column metadata in test_gct
    if (test_gct.col_metadata_df).shape[1] == 0:
        logger.info("No metadata provided for test GCT columns. " +
                    "Using ids as perturbation identifiers.")
        test_gct.col_metadata_df[query_field_name] = test_gct.col_metadata_df.index

    else:
        # If no aggregation fields were provided, use indices
        if len(fields_to_aggregate_in_test_gct_queries) == 0:
            logger.info("No aggregation fields provided for test GCT columns. " +
                        "Using ids as perturbation identifiers.")
            test_gct.col_metadata_df[query_field_name] = test_gct.col_metadata_df.index

        # Otherwise, create a new aggregated field
        else:
            test_gct.col_metadata_df = aggregate_fields(
                test_gct.col_metadata_df,
                fields_to_aggregate_in_test_gct_queries, sep, query_field_name)

    # Check if we have any row metadata in test_gct
    if (test_gct.row_metadata_df).shape[1] == 0:
        logger.info("No metadata provided for test GCT rows. " +
                    "Using ids as perturbation identifiers.")
        test_gct.row_metadata_df[target_field_name] = test_gct.row_metadata_df.index

    else:

        # If no aggregation fields were provided, use indices
        if len(fields_to_aggregate_in_test_gct_targets) == 0:
            logger.info("No aggregation fields provided for test GCT rows. " +
                        "Using ids as perturbation identifiers.")
            test_gct.row_metadata_df[target_field_name] = test_gct.row_metadata_df.index

        # Otherwise, create a new aggregated field
        else:
            test_gct.row_metadata_df = aggregate_fields(
                test_gct.row_metadata_df,
                fields_to_aggregate_in_test_gct_targets, sep, target_field_name)

    # Check if we have any column metadata in bg_gct
    if (bg_gct.col_metadata_df).shape[1] == 0:
        logger.info("No metadata provided for background GCT columns. " +
                    "Using ids as perturbation identifiers.")
        bg_gct.col_metadata_df[target_field_name] = bg_gct.col_metadata_df.index

    else:

        # If no aggregation fields were provided, use indices
        if len(fields_to_aggregate_in_bg_gct) == 0:
            logger.info("No aggregation fields provided for background GCT columns. " +
                        "Using ids as perturbation identifiers.")
            bg_gct.col_metadata_df[target_field_name] = bg_gct.col_metadata_df.index

        # Otherwise, create a new aggregated field
        else:
            bg_gct.col_metadata_df = aggregate_fields(
                bg_gct.col_metadata_df,
                fields_to_aggregate_in_bg_gct, sep, target_field_name)

    # Check if we have any row metadata in bg_gct
    if (bg_gct.row_metadata_df).shape[1] == 0:
        logger.info("No metadata provided for background GCT rows. " +
                    "Using ids as perturbation identifiers.")
        bg_gct.row_metadata_df[target_field_name] = bg_gct.row_metadata_df.index

    else:

        # If no aggregation fields were provided, use indices
        if len(fields_to_aggregate_in_bg_gct) == 0:
            logger.info("No aggregation fields provided for background GCT rows. " +
                        "Using ids as perturbation identifiers.")
            bg_gct.row_metadata_df[target_field_name] = bg_gct.row_metadata_df.index

        # Otherwise, create a new aggregated field
        else:
            bg_gct.row_metadata_df = aggregate_fields(
                bg_gct.row_metadata_df,
                fields_to_aggregate_in_bg_gct, sep, target_field_name)

    return test_gct, bg_gct


def aggregate_fields(df, list_of_fields, separator, agg_field_name):
    """ Join multiple columns of a dataframe into a single column.

    Args:
        df (pandas df)
        list_of_fields (list of strings): columns from df to aggregate
        separator (string): string separator to use in the aggregation
        agg_field_name (string): what to name the new, aggregated column

    Returns:
        df (pandas df): same as input, but one new column appended

    """
    # Check that each field to aggregate is present
    for f in list_of_fields:
        assert f in df.columns, (
            "{} is not present as a metadata field.".format(f))

    sub_df = df[list_of_fields]

    agg_col = []
    for _, my_series in sub_df.iterrows():
        agg_col.append(my_series.astype(str).str.cat(sep=separator))

    df[agg_field_name] = agg_col

    return df


def check_symmetry(test_df, bg_df):
    """ Check whether test and background dfs are symmetric. Assumes matrix
    is symmetric if it is square; this is a soft check! This matters because
    it affects data extraction later (to avoid double-counting).

    Currently, background matrix MUST be square.

    Args:
        test_df (pandas df)
        bg_df (pandas df)

    Returns:
        is_test_df_sym (bool)
        is_bg_df_sym (bool)

    """
    # Check if test_df is square
    if test_df.shape[0] == test_df.shape[1]:
        is_test_df_sym = True
        logger.info("test_df is square, so it will be considered symmetric.")
    else:
        is_test_df_sym = False

    # Check if bg_df is symmetric
    if bg_df.shape[0] == bg_df.shape[1]:
        is_bg_df_sym = True
        logger.info("bg_df is square, so it will be considered symmetric.")
    else:
        is_bg_df_sym = False

    # TODO(lev): should be able to support a non-symmetric background matrix
    # (when we do this, we should also separate bg_gct_field into
    # bg_gct_query_field and bg_gct_target_field)
    assert is_bg_df_sym, "bg_df must be symmetric!"

    return is_test_df_sym, is_bg_df_sym


def compute_connectivities(test_gct, bg_gct, test_gct_query_field,
                           test_gct_target_field, bg_gct_field,
                           connectivity_metric, is_test_df_sym, sep):
    """ Compute all connectivities for a single test_gct and a single bg_gct.

    Args:
        test_gct (GCToo): m rows x n cols, where n is the # of queries, m is the # of targets
        bg_gct (GCToo): M rows x M rows, where M is a superset of m
        test_gct_query_field (string)
        test_gct_target_field (string)
        bg_gct_field (string)
        connectivity_metric (string)
        is_test_df_sym (bool)
        sep (string): separator to use in creating aggregated strings

    Returns:
        conn_gct (GCToo): m rows x n cols, where n is the # of queries, m is the # of targets
        signed_conn_gct (GCToo): m rows x n cols, where n is the # of queries, m is the # of targets

    """
    logger.info("Computing connectivities...")

    # Extract queries from test_df columns and targets from test_df index
    queries = test_gct.col_metadata_df[test_gct_query_field].unique()

    # Extract targets from both test_df and bg_df
    test_targets = test_gct.row_metadata_df[test_gct_target_field].unique()
    bg_targets = bg_gct.row_metadata_df[bg_gct_field].unique()

    # Get intersection of targets
    targets = np.intersect1d(test_targets, bg_targets)
    assert targets.size > 0, (
        "There are no targets in common between the test and bg dfs.\n" +
        "test_targets: {}, bg_targets: {}".format(test_targets, bg_targets))

    logger.info("{} queries, {} targets".format(len(queries), len(targets)))

    # Initialize conn_df and signed_conn_df :: len(targets) x len(queries)
    conn_df = pd.DataFrame(np.nan, index=targets, columns=queries)
    signed_conn_df = conn_df.copy()

    for target in targets:

        # Extract background values
        bg_vals = extract_bg_vals_from_sym(target, bg_gct_field, bg_gct)

        # Make sure bg_vals has at least 1 element before continuing
        # and that bg_vals is not all NaN
        if len(bg_vals) > 0 and not all(np.isnan(bg_vals)):

            for query in queries:
                logger.debug("query: {}, target: {}".format(query, target))

                # Extract test values
                test_vals = extract_test_vals(
                    query, target, test_gct_query_field,
                    test_gct_target_field, test_gct, is_test_df_sym)

                # Make sure test_vals has at least 1 element before continuing
                # and that test_vals is not all NaN
                if len(test_vals) > 0 and not all(np.isnan(test_vals)):

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

    # Aggregate all metadata from test_gct before inserting into output GCT
    col_meta_df_tmp = aggregate_metadata(test_gct.col_metadata_df,
                                         test_gct_query_field, sep)
    row_meta_df_tmp = aggregate_metadata(test_gct.row_metadata_df,
                                         test_gct_target_field, sep)

    # Sort indices of all matrices to ensure that our metadata dfs are in
    # sync with the data dfs
    col_meta_df = col_meta_df_tmp.sort_index(axis=0)
    row_meta_df = row_meta_df_tmp.sort_index(axis=0)
    conn_df_sorted = conn_df.sort_index(axis=0).sort_index(axis=1)
    signed_conn_df_sorted = signed_conn_df.sort_index(axis=0).sort_index(axis=1)

    assert col_meta_df.index.equals(signed_conn_df_sorted.columns), (
        "col_meta_df aggregated entries are not the same as the " +
        "aggregated entries in data_df. col_meta_df.index: {}\n" +
        "signed_conn_df_sorted.columns: {}".format(
            col_meta_df.index, signed_conn_df_sorted.columns))

    assert row_meta_df.index.equals(signed_conn_df_sorted.index), (
        "row_meta_df aggregated entries are not the same as the " +
        "aggregated entries in data_df. row_meta_df.index: {}\n" +
        "signed_conn_df_sorted.index: {}".format(
            row_meta_df.index, signed_conn_df_sorted.index))

    # Create output GCTs
    conn_gct = GCToo.GCToo(data_df=conn_df_sorted,
                           row_metadata_df=row_meta_df,
                           col_metadata_df=col_meta_df)

    signed_conn_gct = GCToo.GCToo(data_df=signed_conn_df_sorted,
                                  row_metadata_df=row_meta_df,
                                  col_metadata_df=col_meta_df)

    return conn_gct, signed_conn_gct


def aggregate_metadata(df, agg_field, sep):
    """ Collapse replicates in a metadata df. agg_field is the column
    metadata field that specifies which rows are replicates of each other.

    Args:
        df (pandas df):
        agg_field (string): column header that contains information about
            which rows should be collapsed with each other; unique version
            of this column becomes the index of out_df
        sep (string): separator to use in creating aggregated strings

    Returns:
        out_df (pandas df): rows (potentially) collapsed

    """
    assert agg_field in df.columns.values, (
        "agg_field {} must be in the columns of df. df.columns.values: {}").format(
            agg_field, df.columns.values)

    # Use lambda so we can pass sep to the function
    out_df = df.groupby(agg_field).agg(lambda x: aggregate_one_series_uniquely(x, sep))

    return out_df


def aggregate_one_series_uniquely(this_series, sep):
    """ String concatenate UNIQUE entries in this_series.

    Args:
        this_series (pandas Series):
        sep (string):

    Returns:
        out_str (string)

    """
    uniq_sorted_str_array = np.sort(this_series.unique()).astype(str)
    out_str = sep.join(uniq_sorted_str_array)

    return out_str


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


def extract_test_vals(query, target, query_field, target_field, test_gct, is_test_df_sym):
    """ Extract values that has query in the columns and target in the rows.

    Args:
        query (string)
        target (string)
        query_field (string): name of metadata field in which to find query
        target_field (string): name of metadata field in which to find target
        test_gct (GCToo)
        is_test_df_sym (bool): only matters if query == target; set to True to
            avoid double-counting in the case of a symmetric matrix

    Returns:
        vals (numpy array)

    """
    assert query in test_gct.col_metadata_df[query_field].values, (
        "query {} is not in the {} metadata of the columns of test_gct.").format(
            query, query_field)

    assert target in test_gct.row_metadata_df[target_field].values, (
        "target {} is not in the {} metadata of the rows of test_gct.").format(
            target, target_field)

    # Extract elements where query is in columns and target is in rows
    target_in_rows_query_in_cols_df = test_gct.data_df.loc[
            test_gct.row_metadata_df[target_field] == target,
            test_gct.col_metadata_df[query_field] == query]

    # If query == target AND the matrix is symmetric, need to take only triu
    # of the extracted values in order to avoid double-counting
    if query == target and is_test_df_sym:
        mask = np.triu(np.ones(target_in_rows_query_in_cols_df.shape), k=1).astype(np.bool)
        vals_with_nans = target_in_rows_query_in_cols_df.where(mask).values.flatten()
        vals = vals_with_nans[~np.isnan(vals_with_nans)]

    else:
        vals = target_in_rows_query_in_cols_df.values.flatten()

    return vals


def extract_bg_vals_from_sym(target, target_field_name, bg_gct):
    """ Extract all values that have some interaction with target.

    Diagonal and lower-right triangle are excluded.

    Args:
        target (string)
        target_field_name (string)
        bg_gct (GCToo)

    Returns:
        vals (numpy array)

    """
    assert target in bg_gct.row_metadata_df[target_field_name].values, (
        "target {} is not in the {} metadata of the rows of bg_gct.".format(
            target, target_field_name))

    assert target in bg_gct.col_metadata_df[target_field_name].values, (
        "target {} is not in the {} metadata of the columns of bg_gct.".format(
            target, target_field_name))

    # Get values in target_field_name in rows
    index_names = bg_gct.row_metadata_df[target_field_name]

    # Get values in target_field_name in columns
    column_names = bg_gct.col_metadata_df[target_field_name]

    # Initialize list of idxs for the values we want to extract
    idxs = []

    # Loop over the rows
    for row_idx in range(bg_gct.data_df.shape[0]):

        # Loop over the columns
        for col_idx in range(bg_gct.data_df.shape[1]):

            # Extract element if it is in the upper-right triangle, excluding
            # the diagonal, and it contains an interaction with target
            if row_idx < col_idx and (index_names[row_idx] == target or column_names[col_idx] == target):
                idxs.append([row_idx, col_idx])

    vals = bg_gct.data_df.values[tuple(np.transpose(idxs))]

    return vals


def extract_bg_vals_from_non_sym(target, target_field_name, bg_gct):
    """ Extract all values that have some interaction with target.

    Just get all values where target is in the rows.

    Args:
        target (string)
        target_field_name (string)
        bg_gct (GCToo)

    Returns:
        vals (numpy array)

    """
    assert target in bg_gct.row_metadata_df[target_field_name].values, (
        "target {} is not in the {} metadata of the rows of bg_gct.".format(
            target, target_field_name))

    vals = bg_gct.data_df.loc[bg_gct.row_metadata_df[target_field_name] == target, :].values.flatten()

    return vals


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