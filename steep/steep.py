import logging
import utils.setup_logger as setup_logger
import argparse
import os
import pandas as pd
import numpy as np
import itertools
from scipy import stats

import utils.psp_utils as utils
import in_out.GCToo as GCToo
import in_out.write_gctoo as write_gctoo

"""
Input is a gct file.

Minimum output is 2 gct files: 1 for similarities and 1 for connectivities.
Could produce more output gct files with alternate metrics related to
connectivity using the "all_conn_outputs" argument.

Similarities are computed pairwise between all samples in in_gct.
Connectivities, however, are aggregated according to "fields." For example,
if fields = ["pert_iname", "cell_id"], perturbations of the same perturbagen in
the same cells (but potentially at different time points and different doses)
are aggregated together in computing connectivity.

"""

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# TODO(lev): move to config file?
OUT_SIM_NAME_SUFFIX = ".steep.similarity.gct"
OUT_CONN_NAME_SUFFIX = ".steep.connectivity.gct"
OUT_PVAL_NAME_SUFFIX = ".steep.pvalue.gct"
OUT_DIR_CONN_NAME_SUFFIX = ".steep.dconnectivity.gct"


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
    parser.add_argument("-fields", nargs="+", default=["pert_iname", "cell_id", "pert_time"],
                        help="list of metadata headers to use for determining unique perturbations")
    parser.add_argument("-out_sim_name", type=str, default=None,
                        help="name of output similarity gct (if None, will use OUT_SIM_NAME_SUFFIX")
    parser.add_argument("-out_conn_name", type=str, default=None,
                        help="name of output connectivity gct (if None, will use OUT_CONN_NAME_SUFFIX")
    parser.add_argument("-out_pval_name", type=str, default=None,
                        help="name of output connectivity p-values gct (if None, will use OUT_PVAL_NAME_SUFFIX")
    parser.add_argument("-out_directed_conn_name", type=str, default=None,
                        help="name of output directed connectivity gct (if None, will use OUT_DIR_CONN_SUFFIX")
    parser.add_argument("-psp_config_path", type=str,
                        default="example_psp.cfg",
                        help="filepath to PSP config file")
    parser.add_argument("-all_conn_outputs", type=str, default=False,
                        action="store_true",
                        help="True will produce more than 1 connectivity gct file")
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                        help="True increases the number of messages reported")

    return parser


def main(args):
    # Get full file path for in_gct, and make sure path exists
    full_in_gct_path = os.path.expanduser(in_gct)
    assert os.path.exists(full_in_gct_path), (
        "The following gct_path cannot be found. full_in_gct_path: {}").format(
            full_in_gct_path)

    # Read gct and config file
    (gct, config_io, config_metadata, config_parameters) = (
        utils.read_gct_and_config_file(in_gct, psp_config))

    # Compute similarity
    similarity_df = compute_similarity_matrix(gct.data_df, method='spearman')

    # Assemble similarity to gct
    similarity_gct = GCToo.GCToo(data_df=similarity_df,
                                 row_metadata_df=gct.col_metadata_df,
                                 col_metadata_df=gct.col_metadata_df)

    # Save similarity gct to output file
    out_sim_name = configure_out_gct_name(full_in_gct_path, OUT_SIM_NAME_SUFFIX, args.out_sim_name)
    full_out_sim_name = os.path.join(out_dir, out_sim_name)
    write_gctoo.write(similarity_gct, full_out_sim_name, filler_null="NA")

    # Compute connectivity
    (conn_df_unpivoted, conn_meta_df) = compute_connectivity(similarity_gct, fields)

    # Subset and pivot out_df_unpivoted to create output dfs
    ks_stat_df = conn_df_unpivoted.pivot("query", "target", "ks_statistic")
    pval_df = conn_df_unpivoted.pivot("query", "target", "p_value")
    ks_stat_directed_df = conn_df_unpivoted.pivot("query", "target", "ks_statistic_directed")

    logger.debug("ks_stat_df.shape: {}".format(ks_stat_df.shape))
    logger.debug("conn_meta_df.shape: {}".format(conn_meta_df.shape))

    # Assemble connectivity to gct
    connectivity_gct = GCToo.GCToo(data_df=ks_stat_df,
                                   row_metadata_df=conn_meta_df,
                                   col_metadata_df=conn_meta_df)

    # Save connectivity gct to output file
    out_conn_name = configure_out_gct_name(full_in_gct_path, OUT_CONN_NAME_SUFFIX, args.out_conn_name)
    full_out_conn_name = os.path.join(out_dir, out_conn_name)
    write_gctoo.write(connectivity_gct, full_out_conn_name, filler_null="NA")

    if args.all_conn_outputs:
        # Assemble p-values to gct
        pval_gct = GCToo.GCToo(data_df=pval_df,
                               row_metadata_df=conn_meta_df,
                               col_metadata_df=conn_meta_df)

        # Save connectivity gct to output file
        out_pval_name = configure_out_gct_name(full_in_gct_path, OUT_PVAL_NAME_SUFFIX, args.out_pval_name)
        full_out_pval_name = os.path.join(out_dir, out_pval_name)
        write_gctoo.write(pval_gct, full_out_pval_name, filler_null="NA")

        # Assemble directed connectivities to gct
        dir_conn_gct = GCToo.GCToo(data_df=ks_stat_directed_df,
                                   row_metadata_df=conn_meta_df,
                                   col_metadata_df=conn_meta_df)

        # Save connectivity gct to output file
        out_dir_conn_name = configure_out_gct_name(full_in_gct_path, OUT_DIR_CONN_NAME_SUFFIX, args.out_directed_conn_name)
        full_out_dir_conn_name = os.path.join(out_dir, out_dir_conn_name)
        write_gctoo.write(dir_conn_gct, full_out_dir_conn_name, filler_null="NA")

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
    # TODO(lev): wtcs?

    return out_df


def compute_connectivity(sim_gct, fields, is_symmetric=True):
    """Computes connectivity of perturbagens with the KS-test statistic.

    Input is a square dataframe of similarity values. Output is a square
    dataframe of connectivity values. The size of the output df will
    be less than or equal to the size of the input df.

    Lower triangle of the similarity matrix is used.

    Connectivity of A to B is determined by a two sample KS-test between
    the distribution of similarities of A to B (all pairwise replicate
    comparisons) and the distribution of similarities of all other
    perturbagens to B.

    Args:
        sim_gct (pandas df): size of data_df = n x n
        fields (list or numpy array of strings): specifies how perturbations are aggregated
        is_symmetric (bool): True indicates that similarity metric is symmetric;
            changes how connectivity is computed

    Returns:
        out_df_unpivoted (pandas df): size = n2 x 5; the columns are query,
            target, ks_statistic, p_value, ks_statistic_directed
        out_meta_df (pandas): row index is perturbation ids, columns are
            metadata headers

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

    # Create copy of sim_gct.data_df on which computations will be performed
    data_df = sim_gct.data_df

    # Aggregate samples according to fields and create a perturbation id
    # for each samples; fields are extracted from row_metadata_df but could be
    # extracted from col_metadata_df just as well
    (pert_ids, out_meta_df) = create_ids_from_metadata(sim_gct.row_metadata_df, "column", fields)

    # Rename ids of data_df to pert_ids
    data_df.index = pd.Index(pert_ids)
    data_df.columns = pd.Index(pert_ids)

    # Set diagonal of similarity matrix to NaN to ignore self-similarities
    np.fill_diagonal(data_df.values, np.nan)

    # if is_symmetric:
    #     # Set upper triangle (incl. diagonal) to NaN
    #     mask = np.tril(np.ones(data_df.shape), k=-1).astype(np.bool)
    #     data_df = data_df.where(mask)
    # else:
    #     # TODO(lev): probably should just symmetrize
    #     err_msg = "Assymetric similarity matrix is not currently supported."
    #     logger.error(err_msg)
    #     raise(Exception(err_msg))

    # Create test and null distributions for KS-test
    (queries, targets, tests, nulls, unique_perts) = create_distributions_for_ks_test(data_df)

    # Initialize output arrays
    ks_stats = np.zeros(np.square(len(unique_perts))) * np.nan
    ks_stats_directed = np.zeros(np.square(len(unique_perts))) * np.nan
    pvals = np.zeros(np.square(len(unique_perts))) * np.nan

    # Perform KS-test between test and null distributions
    count = 0
    logger.info("Computing connectivities...")
    for test, null in zip(tests, nulls):
        try:
            (ks_stat, pval) = stats.ks_2samp(test, null)
        except ValueError as e:
            warning_msg = (
                ("KS-statistic could not be computed for count {}. " +
                 "Error message was the following: {}").format(count, e))
            logger.warning(warning_msg)
            count = count + 1
            continue

        # If median of test distribution is >= median of null
        # distribution, ks_stat_directed = ks_stat; otherwise,
        # ks_stat_directed = ks_stat * -1
        if (np.median(test) - np.median(null)) >= 0:
            ks_stat_directed = ks_stat
        else:
            ks_stat_directed = ks_stat * -1

        # Populate output arrays
        ks_stats[count] = ks_stat
        pvals[count] = pval
        ks_stats_directed[count] = ks_stat_directed

        # Increment counter
        count = count + 1

    # Assemble output arrays into out_df_unpivoted
    out_df_unpivoted = pd.DataFrame.from_dict({
        "query": queries, "target": targets,
        "ks_statistic": ks_stats, "p_value": pvals,
        "ks_statistic_directed": ks_stats_directed})

    # Rearrange columns of out_df_unpivoted
    out_df_unpivoted = out_df_unpivoted[["query", "target", "ks_statistic",
                                         "p_value", "ks_statistic_directed"]]

    return out_df_unpivoted, out_meta_df


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
    # Determine how many unique perturbations there are
    unique_perts = np.unique(in_df.columns.values)
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
        # logger.debug("query: {}, target: {}".format(query, target))

        # Can only proceed if more than 1 target replicate;
        # KS-test cannot be computed otherwise
        num_target_replicates = in_df.columns.isin([target]).sum()

        if num_target_replicates > 1:

            # Extract all elements where target is in the rows
            target_df = in_df.loc[target, :]

            # Extract elements where index and column are both target,
            # set upper triangle of them to NaN, and reinsert into target_df;
            # this is necessary in order to avoid double-counting
            target_equals_query_df = in_df.loc[target, target]
            mask = np.tril(np.ones(target_equals_query_df.shape)).astype(np.bool)
            masked_target_equals_query_df = target_equals_query_df.where(mask)
            # logger.debug("masked_target_equals_query_df.values: {}".format(
            #     masked_target_equals_query_df.values))
            # logger.debug("target_df.loc[target, target].values: {}".format(
            #     target_df.loc[target, target].values))
            target_df.loc[target, target] = masked_target_equals_query_df

            # Elements where column=query is test distribution; remove NaNs
            test = target_df.loc[:, query].values.flatten()
            test = test[~np.isnan(test)]

            # Elements where column!=query is null distribution; remove NaNs
            null = target_df.loc[:, ~target_df.columns.isin([query])].values.flatten()
            null = null[~np.isnan(null)]

            # test must have more than min_distribution_elements elements
            if len(test) < min_distribution_elements:
                warning_msg = (
                    ("For query {} and target {}, there are fewer than {} " +
                    "elements in the test distribution. Connectivity cannot " +
                     "be returned for this combination. test: {}").format(
                        query, target, min_distribution_elements, test))
                logger.warning(warning_msg)
                count = count + 1
                continue

            # null must have more than min_distribution_elements
            if len(null) < min_distribution_elements:
                warning_msg = (
                    ("For query {} and target {}, there are fewer than {} " +
                    "elements in the null distribution. Connectivity cannot " +
                     "be returned for this combination. null: {}").format(
                        query, target, min_distribution_elements, null))
                logger.warning(warning_msg)
                continue

            nulls[count] = null
            tests[count] = test

        # Add increment
        count = count + 1

    return queries, targets, tests, nulls, unique_perts


def symmetrize_if_needed(sim_df):
    """Check if sim_df is symmetric. If not, symmetrize it.

    Args:
        sim_df (pandas df)
    Returns:
        out_df (pandas df): same size as sim_df
    """
    is_sym = np.allclose(sim_df, sim_df.transpose(), equal_nan=True)
    if not is_sym:
        out_df = (sim_df + sim_df.transpose()) / 2
    else:
        out_df = sim_df
    return out_df


def create_ids_from_metadata(metadata_df, dim, fields, sep="_"):
    """Create an id for each row/column using the metadata specified by fields.
    Fields are the metadata fields to use for aggregating perturbagens.

    If dim="row", this fcn will look for fields in the row names and create
    an id for each column. If dim="column", this fcn will look for fields
    in the column names and create an id for each row.

    This fcn also returns subset_df, in which only "fields" remain,
    the index (if dim=column) or columns (if dim=row) are set to ids, and
    duplicate rows (if dim=column) or columns (if dim=row) are removed.

    Args:
        metadata_df (pandas df): size = m x n
        dim (string): choices = {"row", "column"}
        fields (list or numpy array of strings)
        sep (string): separator to use in creating ids

    Returns:
        ids (list of strings): length = m (if dim="col") or n (if dim="row")
        subset_df (pandas df): index (if dim=column) or columns (if dim=row)
            set to ids, duplicate entries removed

    """
    assert (dim == "row" or dim == "column"), (
        "Dim must be either 'row' or 'column'. dim: {}".format(dim))

    # Check that each field exists in metadata_df
    for field in fields:
        if dim == "row":
            assert field in metadata_df.index.values
        elif dim == "column":
            assert field in metadata_df.columns.values

    if dim == "row":
        # Subset df to just the desired rows
        subset_df = metadata_df.loc[fields, :]

        # Convert col to id
        ids = [sep.join(subset_df[col]) for col in subset_df]

        # Set column index of subset_df to ids
        subset_df.columns = pd.Index(ids)

        # Remove duplicate columns
        # (tranpose to use drop_duplicates method and transpose back)
        subset_df = subset_df.transpose().drop_duplicates().transpose()

    elif dim == "column":
        # Subset df to just the desired cols
        subset_df = metadata_df[fields]

        # Convert row to id
        ids = [sep.join(row) for row in subset_df.values]

        # Set index of subset_df to ids
        subset_df.index = pd.Index(ids)

        # Remove duplicate rows
        subset_df = subset_df.drop_duplicates()


    return ids, subset_df


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