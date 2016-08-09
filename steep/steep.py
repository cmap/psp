import sys
import os
import logging
import pandas as pd
import argparse

import utils.setup_logger as setup_logger
import utils.psp_utils as psp_utils
import in_out.GCToo as GCToo
import in_out.parse_gctoo as pg
import in_out.write_gctoo as wg

"""
steep.py

This module computes similarity (Pearson or Spearman correlation). If given 1
gct, steep will compute all pairwise similarities between its columns. If given
2 gcts, steep will compute pairwise similarities between the columns of gct1
and the columns of gct2, but not the similarities within a single gct.

Required arguments are 1 gct, the output directory, and the name for the
output similarity gct.

"""

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("in_gct", type=str, help="path to input gct file")
    parser.add_argument("out_dir", type=str, help="path to output directory")
    parser.add_argument("out_name", type=str, help="what to name the output similarity file")

    # Optional args
    parser.add_argument("-in_gct2", type=str, help="path to second gct file")
    parser.add_argument("-similarity_metric", type=str, default="spearman",
                        choices=["spearman", "pearson"],
                        help="similarity metric to use for comparing columns")
    parser.add_argument("-psp_config_path", type=str,
                        default="psp_production.cfg",
                        help="filepath to PSP config file")
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser

def main(args):

    # Read in the first gct
    (gct1, _, _, _) = psp_utils.read_gct_and_config_file(args.in_gct, args.psp_config_path)

    # If second gct provided, compute similarity between 2 gcts
    if args.in_gct2 is not None:
        logger.info("in_gct2 was provided. Will compute pairwise similarities " +
                    "between the columns of in_gct and in_gct2.")

        # Read in the second gct
        (gct2, _, _, _) = psp_utils.read_gct_and_config_file(args.in_gct2, args.psp_config_path)

        # Compute similarities between gct1 and gct2
        out_df = compute_similarity_bw_two_dfs(gct1.data_df, gct2.data_df, args.similarity_metric)

        # Row metadata is from gct1, column metadata is from gct2
        row_metadata_df = gct1.col_metadata_df
        col_metadata_df = gct2.col_metadata_df

        # Assemble output gct
        out_gct = GCToo.GCToo(out_df, row_metadata_df, col_metadata_df)

    # If only 1 gct provided, compute similarities between the columns of gct1
    else:
        out_df = compute_similarity_within_df(gct1.data_df, args.similarity_metric)

        # Row and column metadata are both from gct1
        metadata_df = gct1.col_metadata_df

        # Assemble output gct
        out_gct = GCToo.GCToo(out_df, metadata_df, metadata_df)

    # Write output gct
    correct_out_name = wg.append_dims_and_file_extension(args.out_name, out_df)
    full_out_name = os.path.join(args.out_dir, correct_out_name)
    wg.write(out_gct, full_out_name, data_null="NaN", metadata_null="NaN", filler_null="NaN")


def compute_similarity_bw_two_dfs(df1, df2, similarity_metric):
    """ Compute similarity between the columns of df1 and the columns of df2.

    The dfs are concated, all pairwise similarities are computed, and then only
    the requested ones (namely between df1 and df2, not within df1  or within
    df2) are returned. This can almost certainly be implemented more
    efficiently, but this method is faster than iterating over columns with a
    for-loop.

    Args:
        df1 (pandas df): size = m x n1
        df2 (pandas df): size = m x n2
        similarity_metric (string): "pearson" or "spearman"

    Returns:
        out_df (pandas df): size = n1 x n2

    """
    # Concatenate the matrices together
    df_concat = pd.concat([df1, df2], axis=1)

    # Keep track of which elements will correspond to the requested output
    df1_cols = range(df1.shape[1])
    df2_cols = range(df1.shape[1], df1.shape[1] + df2.shape[1])

    # Compute similarity
    if similarity_metric is "pearson":
        full_df = df_concat.corr(method="pearson")
    elif similarity_metric is "spearman":
        full_df = df_concat.corr(method="spearman")
    else:
        err_msg = ("similarity metric must be 'pearson' or 'spearman'. " +
                   "similarity_metric: {}").format(similarity_metric)
        raise(Exception(err_msg))

    # Return just the subset of data that was asked for
    out_df = full_df.iloc[df1_cols, df2_cols]

    return out_df


def compute_similarity_within_df(df, similarity_metric):
    """ Compute all pairwise similarities between the columns of df.

    Args:
        df (pandas df): size = m x n
        similarity_metric (string): "pearson" or "spearman"

    Returns:
        out_df (pandas df): size = n x n

    """
    # Compute similarity
    if similarity_metric is "pearson":
        out_df = df.corr(method="pearson")
    elif similarity_metric is "spearman":
        out_df = df.corr(method="spearman")
    else:
        err_msg = ("similarity metric must be 'pearson' or 'spearman'. " +
                   "similarity_metric: {}").format(similarity_metric)
        raise(Exception(err_msg))

    return out_df


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose = args.verbose)

    main(args)