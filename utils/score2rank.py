"""

For each file within in_dir that matches file_wildcard, the score file will
be converted to a rank file. If the input name is INPUT_SEPARATOR_SCORE.gct,
the output name for that file will be INPUT_RANK.gct.

"""

import sys
import argparse
import glob
import os
import numpy as np

import GCToo.GCToo as GCToo
import GCToo.parse_gctoo as pg
import GCToo.write_gctoo as wg


def build_parser():
    """ Build argument parser. """

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--in_dir", "-i", required=True,
                        help="path to input directory")
    parser.add_argument("--out_dir", "-o", required=True,
                        help="path to output directory")
    parser.add_argument("--file_wildcard", "-w", required=True,
                        help=("wildcard to use for finding files within in_dir " +
                              "(make sure to surround in quotes when calling from command line!)"))

    # Optional args
    parser.add_argument("--prefix_separator", "-p", default="_cs_",
                        help=("separator for spltting the file name in order " +
                              "to create a sensible output file name"))
    parser.add_argument("--output_suffix", "-s", default="_RANK.gct",
                        help=("separator for spltting the file name in order " +
                              "to create a sensible output file name"))

    return parser


def main(args):

    # Find files
    full_path_wildcard = args.in_dir + args.file_wildcard
    gct_paths = glob.glob(full_path_wildcard)

    assert len(gct_paths) > 1, "full_path_wildcard: {}".format(full_path_wildcard)

    # Extract prefixes in order to use them later for saving
    prefixes = [(os.path.basename(path)).split(args.prefix_separator)[0] for path in gct_paths]

    for path, prefix in zip(gct_paths, prefixes):
        print "path: {}".format(path)
        print "prefix: {}".format(prefix)

    # Import gcts
    gctoos = [pg.parse(x) for x in gct_paths]

    assert len(gctoos) > 1, "gct_paths: {}".format(gct_paths)

    # Compute & save ranks
    for g, prefix in zip(gctoos, prefixes):

        # Extract data_df
        score_df = g.data_df

        # Must be square
        assert score_df.shape[0] == score_df.shape[1], "Input dataframe must be square."

        # Set diagonal to NaN
        np.fill_diagonal(score_df.values, np.nan)

        # Rank the matrix
        rank_df = score_df.rank(ascending=False)

        # Make a GCToo
        rank_gctoo = GCToo.GCToo(
            data_df=rank_df,
            row_metadata_df=g.row_metadata_df,
            col_metadata_df=g.col_metadata_df)

        # Save the rank_df to file
        out_name = out_dir + prefix + args.output_suffix
        wg.write(rank_gctoo, out_name, filler_null="NaN", data_null="NaN", metadata_null="NaN")

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    main(args)