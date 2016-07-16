import logging
import utils.setup_logger as setup_logger
import argparse
import sys
import os
import pandas as pd
import numpy as np
import itertools
from scipy import stats

import utils.psp_utils as utils
import in_out.GCToo as GCToo
import in_out.write_gctoo as write_gctoo

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

"""
Input is a gct file.

Minimum output is 2 gct files: 1 for similarities and 1 for connectivities.
Could produce 2 gct files with alternate metrics of connectivity --
p-value and signed KS-test statistic -- using the "all_conn_outputs" argument.

Similarities are computed pairwise between all samples in in_gct.
Connectivities, however, are aggregated according to group_fields. For example,
if group_fields = ["pert_iname", "cell_id"], perturbations of the same
perturbagen in the same cell lines (but potentially at different time points
and different doses) are aggregated together when computing connectivity.

"""

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# TODO(lev): move to config file?
OUT_SIM_NAME_SUFFIX = ".steep.similarity.gct"
OUT_CONN_NAME_SUFFIX = ".steep.connectivity.gct"
OUT_PVAL_NAME_SUFFIX = ".steep.pvalue.gct"
OUT_CONN_SIGNED_NAME_SUFFIX = ".steep.connectivity_signed.gct"


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("in_gct", type=str, help="filepath to input gct file")
    parser.add_argument("out_dir", type=str, help="path to output directory")

    # Optional args
    parser.add_argument("-similarity_method", type=str,
                        default="spearman",
                        choices=["spearman", "pearson"],
                        help="method for computing similarity")
    parser.add_argument("-group_fields", nargs="+", default=["pert_iname", "cell_id", "pert_time"],
                        help="list of metadata headers to use for determining unique perturbations")
    parser.add_argument("-out_sim_name", type=str, default=None,
                        help="name of output similarity gct (if None, will use OUT_SIM_NAME_SUFFIX")
    parser.add_argument("-out_conn_name", type=str, default=None,
                        help="name of output connectivity gct (if None, will use OUT_CONN_NAME_SUFFIX")
    parser.add_argument("-out_pval_name", type=str, default=None,
                        help="name of output connectivity p-values gct (if None, will use OUT_PVAL_NAME_SUFFIX")
    parser.add_argument("-out_conn_signed_name", type=str, default=None,
                        help="name of output signed connectivity gct (if None, will use OUT_CONN_SIGNED_NAME_SUFFIX")
    parser.add_argument("-psp_config_path", type=str,
                        default="psp_production.cfg",
                        help="filepath to PSP config file")
    parser.add_argument("-all_conn_outputs", "-aco", action="store_true", default=False,
                        help="whether to increase the # of gct files produced")
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):
    # Get full file path for in_gct, and make sure path exists
    full_in_gct_path = os.path.expanduser(args.in_gct)
    assert os.path.exists(full_in_gct_path), (
        "The following gct_path cannot be found. full_in_gct_path: {}").format(
            full_in_gct_path)

    # Read gct and config file
    (gct, config_io, config_metadata, config_parameters) = (
        utils.read_gct_and_config_file(args.in_gct, args.psp_config_path))

    ### COMPUTE SIMILARITY
    similarity_df = compute_similarity_matrix(gct.data_df, method="spearman")

    # Assemble similarity to gct
    similarity_gct = GCToo.GCToo(data_df=similarity_df,
                                 row_metadata_df=gct.col_metadata_df,
                                 col_metadata_df=gct.col_metadata_df)

    # Save similarity gct to output file
    out_sim_name = configure_out_gct_name(full_in_gct_path, OUT_SIM_NAME_SUFFIX, args.out_sim_name)
    full_out_sim_name = os.path.join(args.out_dir, out_sim_name)
    write_gctoo.write(similarity_gct, full_out_sim_name, filler_null="NA")

    ### COMPUTE CONNECTIVITY
    (conn_df_unpivoted, conn_meta_df) = compute_connectivity(similarity_gct, args.group_fields)

    # Assemble connectivity outputs to gcts
    (connectivity_gct, pval_gct, signed_conn_gct) = assemble_output_conn_gcts(conn_df_unpivoted, conn_meta_df)

    # Save connectivity gct to output file
    out_conn_name = configure_out_gct_name(full_in_gct_path, OUT_CONN_NAME_SUFFIX, args.out_conn_name)
    full_out_conn_name = os.path.join(args.out_dir, out_conn_name)
    write_gctoo.write(connectivity_gct, full_out_conn_name, filler_null="NA")

    if args.all_conn_outputs:
        # Save p-values gct to output file
        out_pval_name = configure_out_gct_name(full_in_gct_path, OUT_PVAL_NAME_SUFFIX, args.out_pval_name)
        full_out_pval_name = os.path.join(args.out_dir, out_pval_name)
        write_gctoo.write(pval_gct, full_out_pval_name, filler_null="NA")

        # Save signed connectivity gct to output file
        out_dir_conn_name = configure_out_gct_name(full_in_gct_path, OUT_CONN_SIGNED_NAME_SUFFIX, args.out_conn_signed_name)
        full_out_dir_conn_name = os.path.join(args.out_dir, out_dir_conn_name)
        write_gctoo.write(signed_conn_gct, full_out_dir_conn_name, filler_null="NA")

    return similarity_gct, connectivity_gct


def compute_similarity_matrix(data_df, method):
    """Compute all pairwise similarities between columns of input df.

    N.B. Row and column ids are in the same order as the input.

    Args:
        data_df (pandas df): size = m x n
        metric (string): similarity metric to use; choices = {"spearman"}

    Returns:
        out_df (pandas df): size = n x n
    """
    logger.info("Computing similarity with method = {}.".format(method))

    if method.lower() == "pearson":
        out_df = data_df.corr(method="pearson") # automatically deals with missing data
    elif method.lower() == "spearman":
        out_df = data_df.corr(method="spearman") # automatically deals with missing data
    # TODO(lev): wtcs

    return out_df


def compute_connectivity(sim_gct, group_fields, is_symmetric=True):
    """Computes connectivity of perturbagens with the KS-test statistic.

    Input is a square dataframe of similarity values. Output is a square
    dataframe of connectivity values. The size of the output df will
    be less than or equal to the size of the input df.

    Assumed that similarity matrix is symmetric, so only its lower triangle
    is used (upper triangle discarded).

    Connectivity of A to B is determined by a two sample KS-test between
    the distribution of similarities of A to B (all pairwise replicate
    comparisons) and the distribution of similarities of all other
    perturbations to B.

    N.B. query and target columns in out_df_unpivoted are sorted case-insensitively.

    Args:
        sim_gct (pandas df): size of data_df = n x n
        group_fields (list or numpy array of strings): specifies how perturbations are aggregated
        is_symmetric (bool): True indicates that similarity metric is symmetric;
            changes how connectivity is computed

    Returns:
        out_df_unpivoted (pandas df): size = m x 5; m < n; the columns are
            query, target, ks_statistic, p_value, ks_statistic_signed
        conn_meta_df (pandas df): metadata df to be used for building
            connectivity gct; row index is group_ids, column headers are group_fields

    """
    # sim_gct.data_df should be square, have the same index and columns,
    # and should have no missing values
    assert sim_gct.data_df.shape[0] == sim_gct.data_df.shape[1], (
        "sim_gct.data_df must be square. sim_gct.data_df.shape: {}".format(
            sim_gct.data_df.shape))
    assert np.array_equal(sim_gct.data_df.columns.values, sim_gct.data_df.index.values), (
        ("sim_gct.data_df.columns.values and sim_gct.data_df.index.values must be the same." +
         "\nsim_gct.data_df.columns.values:\n{}\nsim_gct.data_df.index.values:\n{}".format(
             sim_gct.data_df.columns.values, sim_gct.data_df.index.values)))
    assert any(sim_gct.data_df.notnull()), "sim_gct.data_df must have no missing values."

    # Create a pert_id for each perturbation using the values in group_fields;
    # fields are extracted from row_metadata_df but could be extracted from
    # col_metadata_df just as well
    (group_ids, meta_df_with_group_ids) = create_group_ids(sim_gct.row_metadata_df, "col", group_fields)

    # Set group_ids to the index of meta_df_with_group_ids and drop_duplicates
    # to produce the appropriate metadata_df for the connectivity gct
    conn_meta_df = meta_df_with_group_ids.set_index("group_id").drop_duplicates()

    # create_distributions_for_ks_test works on a df where the indices are
    # group_ids; this means the indices might not be unique
    data_df_for_creating_distributions = sim_gct.data_df.copy()
    data_df_for_creating_distributions.index = pd.Index(group_ids)
    data_df_for_creating_distributions.columns = pd.Index(group_ids)

    # Set diagonal of similarity matrix to NaN to ignore self-similarities
    np.fill_diagonal(data_df_for_creating_distributions.values, np.nan)

    # Create test and null distributions for KS-test
    logger.info("Creating distributions for KS tests...")
    (queries, targets, tests, nulls, unique_perts) = (
        create_distributions_for_ks_test(data_df_for_creating_distributions))

    # Initialize output arrays
    ks_stats = np.zeros(np.square(len(unique_perts))) * np.nan
    ks_stats_signed = np.zeros(np.square(len(unique_perts))) * np.nan
    pvals = np.zeros(np.square(len(unique_perts))) * np.nan

    # Perform KS-test between test and null distributions
    count = 0
    logger.info("Computing connectivities...")
    for test, null, query, target in zip(tests, nulls, queries, targets):
        # Compute KS-test only if both distributions are not empty
        if len(test)!=0 and len(null)!=0:

            try:
                (ks_stat, pval) = stats.ks_2samp(test, null)
            except ValueError as e:
                warning_msg = (
                    ("KS-statistic could not be computed for query={}, target={}. " +
                     "Error message was the following: {}\ntest:\n{}\nnull:\n{}.").format(
                        query, target, e, test, null))
                logger.warning(warning_msg)
                count = count + 1
                continue

            # If median of test distribution is >= median of null
            # distribution, ks_stat_signed = ks_stat; otherwise,
            # ks_stat_signed = ks_stat * -1
            if (np.median(test) - np.median(null)) >= 0:
                ks_stat_signed = ks_stat
            else:
                ks_stat_signed = ks_stat * -1

            # Populate output arrays
            ks_stats[count] = ks_stat
            pvals[count] = pval
            ks_stats_signed[count] = ks_stat_signed

            # See specific output for debugging purposes
            if query=="DMSO_PC3_24" and target=="DMSO_PC3_24":
                logger.debug("ks_stat: {}, pval: {}".format(ks_stat, pval))

        # Increment counter
        count = count + 1

    # Assemble output arrays into out_df_unpivoted
    out_df_unpivoted = pd.DataFrame.from_dict({
        "query": queries, "target": targets,
        "ks_statistic": ks_stats, "p_value": pvals,
        "ks_statistic_signed": ks_stats_signed})

    # Rearrange columns of out_df_unpivoted
    out_df_unpivoted = out_df_unpivoted[["query", "target", "ks_statistic",
                                         "p_value", "ks_statistic_signed"]]

    return out_df_unpivoted, conn_meta_df


def create_distributions_for_ks_test(in_df, min_distribution_elements=2):
    """Use index and column names to extract test and null
    distributions for each unique perturbation. The number of unique
    perturbations is equal to n, or len(unique_perts).

    Distributions are extracted in the following process: first, all rows
    where target is the index name are extracted. Then, half of the elements
    where index name and column name are target are set to NaN; this makes
    sure that those elements are not extracted twice. Next, the elements where
    the column is query are extracted to become the test distribution. The
    remaining elements (i.e. where column is NOT query) are extracted to
    become the null distribution.

    Args:
        in_df (pandas df): size = M x N
        min_distribution_elements (int): minimum number of elements that must
            be present in all null and test distributions

    Returns:
        queries (numpy array of strings): length = n
        targets (numpy array of strings): length = n
        tests (numpy array of arrays): length = n
        nulls (numpy array of arrays): length = n
        unique_perts (numpy array of strings): length = n

    N.B. n is the number of unique elements in N, so n < N.

    """
    # Turn off warnings about chained assignment (comes up later when
    # modifying target_df before extracting the test and null distributions)
    pd.options.mode.chained_assignment = None

    # Determine unique perturbations
    unique_perts_case_sensitive_sort = np.unique(in_df.columns.values)

    # Re-sort unique_perts in order for it to be sorted case-insensitively
    unique_perts = np.array(sorted(list(unique_perts_case_sensitive_sort), key=str.lower))

    # Determine number of unique perts
    num_unique_perts = len(unique_perts)
    logger.info("Number of unique perturbations is {}.".format(num_unique_perts))

    # Initialize output arrays
    queries = [[] for x in range(np.square(num_unique_perts))]
    targets = [[] for x in range(np.square(num_unique_perts))]
    tests = [[] for x in range(np.square(num_unique_perts))]
    nulls = [[] for x in range(np.square(num_unique_perts))]

    # Iterate over all query, target combinations
    count = 0
    for query, target in itertools.product(unique_perts, repeat=2):
        queries[count] = query
        targets[count] = target

        # Can only proceed if more than 1 target replicate;
        # KS-test cannot be computed otherwise
        num_target_replicates = in_df.columns.isin([target]).sum()

        if num_target_replicates > 1:

            # Extract all elements where target is in the rows
            # N.B. Taking a copy here ensures that in_df is not modified when
            # we modify target_df later
            target_df = in_df.loc[target, :].copy()

            # Extract elements where index and column are both target,
            # set upper triangle of them to NaN, and reinsert into target_df;
            # this is necessary in order to avoid double-counting
            target_equals_query_df = in_df.loc[target, target].copy()
            mask = np.tril(np.ones(target_equals_query_df.shape)).astype(np.bool)
            masked_target_equals_query_df = target_equals_query_df.where(mask)
            target_df.loc[target, target] = masked_target_equals_query_df

            # Elements where column=query is test distribution; remove NaNs
            test = target_df.loc[:, query].values.flatten()
            test = test[~np.isnan(test)]

            # TODO(lev): need to change here to modify null distribution

            # Elements where column!=query is null distribution; remove NaNs
            null = target_df.loc[:, ~target_df.columns.isin([query])].values.flatten()
            null = null[~np.isnan(null)]

            # test must have more than min_distribution_elements elements
            if len(test) < min_distribution_elements:
                warning_msg = (
                    ("For query {} and target {}, there are fewer than {} " +
                    "elements in the test distribution. Connectivity cannot " +
                     "be returned for this combination, so test and null will "
                     "both be set to []. test: {}").format(
                        query, target, min_distribution_elements, test))
                logger.warning(warning_msg)
                count = count + 1
                continue

            # null must have more than min_distribution_elements
            if len(null) < min_distribution_elements:
                warning_msg = (
                    ("For query {} and target {}, there are fewer than {} " +
                    "elements in the null distribution. Connectivity cannot " +
                     "be returned for this combination, so test and null will "
                     "both be set to []. null: {}").format(
                        query, target, min_distribution_elements, null))
                logger.warning(warning_msg)
                count = count + 1
                continue

            # Null and test are only assigned if there are enough elements in
            # both the null and test distributions
            nulls[count] = null
            tests[count] = test

            # See specific output for debugging purposes
            if query=="DMSO_PC3_24" and target=="DMSO_PC3_24":
                logger.debug("num of test elements: {}; num of null elements: {}".format(
                    len(test), len(null)))
                logger.debug("\ntarget_equals_query_df:\n{}".format(target_equals_query_df))
                logger.debug("\nmasked_target_equals_query_df:\n{}".format(masked_target_equals_query_df))
                logger.debug("\ntest:\n{}\nnull:\n{}".format(test, null))

        # Add increment
        count = count + 1

    return queries, targets, tests, nulls, unique_perts


def create_group_ids(metadata_df, new_dim, group_fields, sep="_"):
    """Create a new field called group_id by joining the strings in
    other specified metadata fields.

    If dim="row", a group_id row will be added, and group_fields must exist
    in the row headers. If dim="col", a group_id column will be added, and
    group_fields must exist in the column headers.

    This fcn also returns subset_df, in which only group_fields and group_id
    remain as the row (if dim="row") or column (if dim="col") headers.

    For example:

        meta_df = pd.DataFrame.from_dict(
            {"cell":["a", "b", "c"], "dose":["10", "10", "1"],
             "time":["24", "24", "24"]})
        group_fields = ["cell", "dose"]
        (group_ids, subset_df) = create_group_ids(meta_df, new_dim="col",group_fields=group_fields)

        # group_ids = ["a_10", "b_10", "c_1"]
        # subset_df = pd.DataFrame.from_dict(
            {"cell":["a", "b", "c"], "dose":["10", "10", "1"], "group_id":["a_10", "b_10", "c_1"]})

    Args:
        metadata_df (pandas df): size = m x n
        new_dim (string): choices = {"row", "col"}
        group_fields (list or numpy array of strings)
        sep (string): separator to use in creating group_ids

    Returns:
        group_ids (list of strings): length = m (if dim="col") or n (if dim="row")
        subset_df (pandas df): if dim="row", size = (len(group_fields) + 1) x
            n; if dim="col", size = m x (len(group_fields) + 1)

    """
    assert (new_dim == "row" or new_dim == "col"), (
        "Dim must be either 'row' or 'col'. dim: {}".format(new_dim))

    # Check that each group_field exists in metadata_df
    for field in group_fields:
        if new_dim == "row":
            assert field in metadata_df.index.values, (
                "{} is not present in the metadata row headers.".format(field))
        elif new_dim == "col":
            assert field in metadata_df.columns.values, (
                "{} is not present in the metadata column headers.".format(field))

    if new_dim == "row":
        # Subset df to just the desired rows
        subset_df = metadata_df.loc[group_fields, :]

        # Convert col to group_id
        group_ids = [sep.join(subset_df[col]) for col in subset_df]

        # Add group_id as a new row
        subset_df.loc["group_id", :] = group_ids


    elif new_dim == "col":
        # Subset df to just the desired cols
        subset_df = metadata_df[group_fields]

        # Convert row to group_id
        group_ids = [sep.join(row) for row in subset_df.values]

        # Add group_id as a new col
        subset_df.loc[:, "group_id"] = group_ids

    return group_ids, subset_df


def assemble_output_conn_gcts(conn_df_unpivoted, conn_meta_df):
    """Assembles numerical connectivity values and metadata into output gcts.

    Args:
        conn_df_unpivoted (pandas df): size = m x 5; the columns are
            query, target, ks_statistic, p_value, ks_statistic_signed
        conn_meta_df (pandas df): metadata df to be used for building
            connectivity gcts; row index is group_ids, column headers are group_fields
    Returns:
        connectivity_gct (GCToo object)
        pval_gct (GCToo object)
        signed_conn_gct (GCToo object)
    """
    # Subset and pivot conn_df_unpivoted to create output dfs
    bad_sorted_ks_stat_df = conn_df_unpivoted.pivot("query", "target", "ks_statistic")
    bad_sorted_pval_df = conn_df_unpivoted.pivot("query", "target", "p_value")
    bad_sorted_ks_stat_signed_df = conn_df_unpivoted.pivot("query", "target", "ks_statistic_signed")

    # Pivot sorts the index, so we have to re-sort case-insensitively
    ks_stat_df = bad_sorted_ks_stat_df.reindex(sorted(bad_sorted_ks_stat_df.index, key=lambda x: x.lower()))
    pval_df = bad_sorted_pval_df.reindex(sorted(bad_sorted_pval_df.index, key=lambda x: x.lower()))
    ks_stat_signed_df = bad_sorted_ks_stat_signed_df.reindex(sorted(bad_sorted_ks_stat_signed_df.index, key=lambda x: x.lower()))

    # Pivot also sorts the columns, so we have to re-sort them case-insensitively as well
    ks_stat_df = ks_stat_df[sorted(ks_stat_df.columns, key=lambda x: x.lower())]
    pval_df = pval_df[sorted(pval_df.columns, key=lambda x: x.lower())]
    ks_stat_signed_df = ks_stat_signed_df[sorted(ks_stat_signed_df.columns, key=lambda x: x.lower())]

    # Make sure that conn_meta_df is also sorted case-insensitively
    sorted_conn_meta_df = conn_meta_df.reindex(sorted(conn_meta_df.index, key=lambda x: x.lower()))

    # Assemble connectivity to gct
    connectivity_gct = GCToo.GCToo(data_df=ks_stat_df,
                                   row_metadata_df=sorted_conn_meta_df,
                                   col_metadata_df=sorted_conn_meta_df)

    # Assemble p-values to gct
    pval_gct = GCToo.GCToo(data_df=pval_df,
                           row_metadata_df=sorted_conn_meta_df,
                           col_metadata_df=sorted_conn_meta_df)

    # Assemble signed connectivities to gct
    signed_conn_gct = GCToo.GCToo(data_df=ks_stat_signed_df,
                                  row_metadata_df=sorted_conn_meta_df,
                                  col_metadata_df=sorted_conn_meta_df)

    return connectivity_gct, pval_gct, signed_conn_gct


def configure_out_gct_name(in_gct, suffix, out_name_from_args):
    """If out_name_from_args is None, append suffix to the input
    gct name. Must end in '.gct'.

    Args:
        in_gct (string)
        suffix (string)
        conn_suffix (string)
        out_name_from_args (string)

    Returns:
        out_name (string)
    """
    input_basename = os.path.basename(in_gct)
    if out_name_from_args is None:
        out_name = input_basename + suffix
    else:
        out_name = out_name_from_args

    assert os.path.splitext(out_name)[1] == ".gct", (
        "The output gct name must end with '.gct'. out_name: {}".format(
            out_name))

    return out_name



if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args: {}".format(args))

    main(args)