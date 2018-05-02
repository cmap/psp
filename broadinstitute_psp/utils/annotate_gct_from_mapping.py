"""
annotate_gct_from_mapping.py

Add annotations to a gct file using a mapping tsv file that looks like this:

from_field  to_field1   to_field2
pert1       moa1        target1
pert2       moa2        target2

This script will add "to_field1" and "to_field2" as metadata fields and populate
those fields using the given mapping. "from_field" in the tsv will be mapped to
the argument "gct_from_field", which must already be in the gct file. Note that
the "from_field" header in the mapping tsv file is not used.

"""

import logging
import sys
import argparse
import pandas as pd

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.write_gct as wg

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--path_to_gct", "-i", required=True,
                        help="path to gct being annotated")
    parser.add_argument("--path_to_mapping_tsv", "-m", required=True,
                        help="path to mapping tsv")

    # Optional args
    parser.add_argument("--out_name", "-o", default="annotated.gct",
                        help="what to name the annotated output gct")
    parser.add_argument("--row_and_or_col", "-rc", choices=["row", "col", "both"],
                        default="both",
                        help=("whether to apply the mapping to only the row metadata, " +
                              "only the col metadata, or both"))
    parser.add_argument("--gct_from_field", "-f", default=None,
                        help=("field already in the gct file; if None, ids " +
                              "will be used"))
    parser.add_argument("--missing_entry", "-me", default="NA",
                        help="what to put for an entry if it doesn't have a mapping")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):

    # Parse gct file
    gct = parse.parse(args.path_to_gct)

    # Parse mapping tsv file
    mapping = pd.read_csv(args.path_to_mapping_tsv, sep="\t", index_col=0)

    # Make sure the ids from the mapping file are unique
    duplicated_bool_array = mapping.index.duplicated()
    assert sum(duplicated_bool_array) == 0, (
        "ids in mapping file must be unique. duplicated ids in mapping:\n{}".format(
        mapping.index[duplicated_bool_array]))

    for col in mapping.columns:

        if args.row_and_or_col == "both":
            annotate_meta_df(gct.row_metadata_df, mapping.loc[:, col], args.gct_from_field, args.missing_entry)
            annotate_meta_df(gct.col_metadata_df, mapping.loc[:, col], args.gct_from_field, args.missing_entry)
        elif args.row_and_or_col == "row":
            annotate_meta_df(gct.row_metadata_df, mapping.loc[:, col], args.gct_from_field, args.missing_entry)
        elif args.row_and_or_col == "col":
            annotate_meta_df(gct.col_metadata_df, mapping.loc[:, col], args.gct_from_field, args.missing_entry)

    wg.write(gct, args.out_name, filler_null="NA", data_null="NaN", metadata_null="NA")


def annotate_meta_df(meta_df, to_entries, gct_from_field, missing_entry):
    """ meta_df is modified in-place.

    Args:
        meta_df (pandas df)
        to_entries (pandas series)
        gct_from_field (string): name of field in gct that has the source entries
        missing_entry (string)

    Returns:
        None

    """
    # If no metadata field is given, use index values
    if gct_from_field is None:
        entries_from_gct = meta_df.index.values

    else:
        assert gct_from_field in meta_df.columns, ("gct_from_field must be a metadata header. " +
            "meta_df.columns.values: {}, gct_from_field: {}".format(
            meta_df.columns.values, gct_from_field))

        entries_from_gct = meta_df.loc[:, gct_from_field]

    # Create mapped_metadata as an array
    mapped_metadata = [to_entries.loc[entry_in_gct] if entry_in_gct in to_entries.index else missing_entry for entry_in_gct in entries_from_gct]

    # Insert destination_metadata into row and col metadata
    meta_df[to_entries.name] = mapped_metadata

    return None


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)
