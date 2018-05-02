"""
sim2dist.py

Convert a similarity matrix to a distance matrix.

d = 1 - s

d = distance
s = similarity

"""

import logging
import argparse
import sys

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.write_gct as wg

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--in_gct_path", "-i", required=True,
                        help="path to input gct")
    parser.add_argument("--out_name", "-o", default="out_dist.gct",
                        help="what to name the output gct")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="Whether to print a bunch of output.")

    return parser


def main(args):

    # Import data
    in_gct = parse.parse(args.in_gct_path)

    # Compute distance df
    dist_df = 1 - in_gct.data_df

    # Create distance gct
    dist_gct = GCToo.GCToo(dist_df, in_gct.row_metadata_df, in_gct.col_metadata_df)

    # Write dist_gct to file
    wg.write(dist_gct, args.out_name, filler_null="NA")


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)
