"""
separate_gct.py

Separates a gct into several gcts.

"""
import logging
import sys
import os
import argparse

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.subset_gctoo as sg
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.write_gct as wg
import cmapPy.pandasGEXpress.write_gctx as wgx

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--in_gct_path", "-i", required=True,
                        help="file path to input gct")
    parser.add_argument("--separate_field", "-sf", required=True,
                        help="field along which to separate the gct")
    parser.add_argument("--row_or_col", "-rc", choices=["row", "col"], required=True,
                        help="whether separate_field is in the row or column metadata")

    # Optional args
    parser.add_argument("--out_dir", "-od", default=".", help="directory in which to save output gcts")
    parser.add_argument("--out_name_prefix", "-op", default="", help="prefix for output file names")
    parser.add_argument("--out_name_suffix", "-os", default=".gct", help="suffix for output file names")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):
    """ The main method. """

    # Import gct
    in_gct = parse.parse(args.in_gct_path)

    # Create the separated gcts
    (out_gcts, out_gct_prefixes) = separate(in_gct, args.separate_field, args.row_or_col)

    # Save the returned gcts
    for gct, name in zip(out_gcts, out_gct_prefixes):
        full_out_name = os.path.join(args.out_dir, args.out_name_prefix + str(name) + args.out_name_suffix)

        # Write to GCT or GCTX depending on extension
        if str.lower(os.path.splitext(full_out_name)[1]) == ".gct":
            wg.write(gct, full_out_name, data_null="NaN", metadata_null="NA", filler_null="NA")
        elif str.lower(os.path.splitext(full_out_name)[1]) == ".gctx":
            wgx.write(gct, full_out_name)
	else:
	    raise(Exception("out_name_suffix must end in either .gct or .gctx. out_name_suffix: {}".format(
	         (args.out_name_suffix))))


def separate(in_gct, separate_field, row_or_col):
    """ Create a new GCT object for each unique value in separate_field.

    Args:
        in_gct (GCToo object)
        separate_field (string)
        row_or_col (string)

    Returns:
        gcts (list of GCToo objects)
        unique_values_in_field (list of strings)

    """
    if row_or_col == "row":
        assert separate_field in in_gct.row_metadata_df.columns, (
            ("separate_field must be in in_gct.row_metadata_df.columns. " +
             "separate_field: {}, in_gct.row_metadata_df.columns: {}").format(
                separate_field, in_gct.row_metadata_df.columns.values))

        unique_values_in_field = list(in_gct.row_metadata_df.loc[:, separate_field].unique())

        gcts = []
        for val in unique_values_in_field:
            bool_array = in_gct.row_metadata_df.loc[:, separate_field].values == val

            new_gct = sg.subset_gctoo(in_gct, row_bool=bool_array)
            gcts.append(new_gct)

    elif row_or_col == "col":
        assert separate_field in in_gct.col_metadata_df.columns, (
            ("separate_field must be in in_gct.col_metadata_df.columns. " +
             "separate_field: {}, in_gct.col_metadata_df.columns: {}").format(
                separate_field, in_gct.col_metadata_df.columns.values))

        unique_values_in_field = list(in_gct.col_metadata_df.loc[:, separate_field].unique())

        gcts = []
        for val in unique_values_in_field:
            bool_array = in_gct.col_metadata_df.loc[:, separate_field].values == val
            new_gct = sg.subset_gctoo(in_gct, col_bool=bool_array)
            gcts.append(new_gct)

    else:
        raise(Exception("row or col must be 'row' or 'col'."))

    # Make sure each gct is associated with a value from separate_field
    assert len(gcts) == len(unique_values_in_field), (
        "len(gcts): {}, len(unique_values_in_field): {}".format(
            len(gcts), len(unique_values_in_field)))

    return gcts, unique_values_in_field


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)
