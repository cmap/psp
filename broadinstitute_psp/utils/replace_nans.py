"""
Replace NaNs in the data df of a GCT with either...

0) zero,
1) the probe median, or
2) the probe mean.
"""

import logging
import argparse
import os
import sys

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.subset_gctoo as sg
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.write_gct as wg

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                         help="Whether to print a bunch of output.")
    parser.add_argument("in_gct_path", type=str,
                        help="path to input gct")
    parser.add_argument("out_name", type=str,
                        help="what to name the output gct")
    parser.add_argument("-replace_with", "-rw", choices=["zero", "median", "mean"],
                        help="what to replace NaN with", default="mean")
    
    return parser


def main(args):

    # Import data
    assert os.path.exists(args.in_gct_path), (
        "in_gct_path could not be found: {}").format(args.in_gct_path)
    in_gct = parse.parse(args.in_gct_path)

    # First, check if any rows are all NaN; if so, remove them
    dropped_df = in_gct.data_df.dropna(how="all")
    bools_of_remaining = in_gct.data_df.index.isin(dropped_df.index.values)
    in_gct = sg.subset_gctoo(in_gct, row_bool=bools_of_remaining)

    if args.replace_with == "zero":
        in_gct.data_df.fillna(0, inplace=True)

    elif args.replace_with == "median":
        probe_medians = in_gct.data_df.median(axis=1)

        for row_idx, row in enumerate(in_gct.data_df.values):
            this_row = in_gct.data_df.iloc[row_idx, :]
            this_row[this_row.isnull()] = probe_medians[row_idx]
            in_gct.data_df.iloc[row_idx, :] = this_row

    elif args.replace_with == "mean":
        probe_means = in_gct.data_df.mean(axis=1)

        for row_idx, row in enumerate(in_gct.data_df.values):
            this_row = in_gct.data_df.iloc[row_idx, :]
            this_row[this_row.isnull()] = probe_means[row_idx]
            in_gct.data_df.iloc[row_idx, :] = this_row

    wg.write(in_gct, args.out_name, filler_null="NA")

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)
