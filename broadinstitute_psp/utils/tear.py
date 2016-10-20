""" Z-scores data.

Converts level 3 to level 4 data. Required input is a filepath to a gct file and
a filepath to an output file. Output is writing a processed gct file.

N.B. This script requires a configuration file. You can specify the location
of this config file with the optional argument -psp_config_path. Otherwise,
it will look for the example configuration file (example_psp.cfg) in the
current directory.

TODO(lev) --> Example usage:

"""

import argparse
import logging
import sys

import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.utils.psp_utils as psp_utils
import broadinstitute_psp.dry.dry as dry

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

ZSCORE_PROV_CODE_ENTRY = "ZSC"
PROV_CODE_FIELD = "provenance"

# TODO(lev): fix me up


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("gct_path", type=str, help="filepath to gct file")
    parser.add_argument("out_path", type=str, help="path to save directory")

    # Optional args
    parser.add_argument("-out_name", type=str, default=None,
                        help="name of output gct (if None, will use <INPUT_GCT>.tear.processed.gct")
    parser.add_argument("-psp_config_path", type=str,
                        default="~/psp_production.cfg",
                        help="filepath to PSP config file")
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                        help="increase the number of messages reported")

    return parser


def main(args):
    # Read gct and config file
    (in_gct, config_io, config_metadata, _) = (
        psp_utils.read_gct_and_config_file(args.gct_path, args.psp_config_path))

    # Extract provenance code
    prov_code = psp_utils.extract_prov_code(
        in_gct.col_metadata_df, config_metadata["prov_code_field"],
        config_metadata["prov_code_delimiter"])

    # Robust z-score
    zscored_df = robust_zscore(in_gct.data_df)
    (zscored_gct, updated_prov_code) = dry.update_metadata_and_prov_code(
        zscored_df, in_gct.row_metadata_df, in_gct.col_metadata_df,
        ZSCORE_PROV_CODE_ENTRY, prov_code)

    # Configure output name
    (out_gct_name, _) = dry.configure_out_names(args.gct_path, args.out_name, None)

    # Reinsert provenance code
    out_gct = dry.insert_offsets_and_prov_code(
        zscored_gct, None, None, updated_prov_code,
        config_metadata["prov_code_field"], config_metadata["prov_code_delimiter"])

    # Write output gct
    dry.write_output_gct(out_gct, args.out_path, out_gct_name,
                         config_io["data_null"], config_io["filler_null"])
    return out_gct


def robust_zscore(data_df):
    """Normalize data using robust z-score.

    Formula: x' = (x - row_median) / (row_mad * 1.4826)

    Args:
        data_df (pandas df)
    Returns:
        out_df (pandas df): normalized
    """
    row_medians = data_df.median(axis=1)
    row_mads = (data_df.sub(row_medians, axis=0)).abs().median(axis=1)
    out_df = data_df.sub(row_medians, axis=0).divide(
        1.4826 * row_mads, axis='index')
    return out_df


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args: {}".format(args))

    main(args)