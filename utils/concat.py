"""
concat.py

Concatenate gct files together.

1) get_file_list
2) concat_row_meta
3) concat_col_meta
4) concat_data

"""


import argparse
import os
import sys
import glob
import logging
import utils.setup_logger as setup_logger
import pandas as pd
import numpy as np
#TODO(lev): erase dependencies on numpy

import in_out.parse_gctoo as parse_gctoo
import in_out.GCToo as GCToo
import in_out.write_gctoo as write_gctoo

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("file_wildcard", type=str,
        help=("wildcard specifying where files should be found " +
              "(surround in quotes if calling from command line)"))
    parser.add_argument("full_out_name", type=str,
        help="what to name the output file (full path)")

    parser.add_argument("-fields_to_remove", type=str, nargs="*",
        default=["pr_probe_suitability_manual", "pr_probe_normalization_group"],
        help="fields to remove from the row metadata headers before concatenating")
    parser.add_argument("-reset_sample_ids", "-rsi", action="store_true", default=False,
        help="whether to reset sample ids (use this flag if sample ids are not unique)")
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
        help="whether to print a bunch of output")

    return parser


def main(args):
    # Find files
    files = get_file_list(args.file_wildcard)

    # Parse each file and append to a list
    gctoos = []
    for f in files:
        gctoos.append(parse_gctoo.parse(f))

    # Create concatenated gctoo object
    out_gctoo = concat_cols(gctoos, args.fields_to_remove, args.reset_sample_ids)

    # Write out_gctoo to file
    logger.info("Write to file...")
    write_gctoo.write(out_gctoo, args.full_out_name, filler_null="NA", data_null="NA")

def get_file_list(wildcard):
    """Search for files to be concatenated.

    Currently very basic, but could expand to be more sophisticated.

    Args:
        wildcard (regular expression string)
    Returns:
        files (list of full file paths)
    """
    files = glob.glob(os.path.expanduser(wildcard))
    return files


def concat_cols(gctoos, fields_to_remove, reset_sample_ids):

    # Separate each gctoo into its component dfs
    row_meta_dfs = []
    col_meta_dfs = []
    data_dfs = []
    for g in gctoos:
        row_meta_dfs.append(g.row_metadata_df)
        col_meta_dfs.append(g.col_metadata_df)
        data_dfs.append(g.data_df)

    # Concatenate row metadata
    all_row_metadata_df = concat_row_meta(row_meta_dfs, fields_to_remove)

    # Concatenate col metadata
    # all_col_metadata_df = concat_col_meta(col_meta_dfs, reset_sample_ids)

    # TODO(lev): keep going here

    # Concatenate data

    # Concatenate the col_metadata_dfs
    all_col_metadata_df = pd.concat([x.col_metadata_df for x in gctoos], axis=0)
    logger.info("all_col_metadata_df.shape: {}".format(all_col_metadata_df.shape))

    # Concatenate the data_dfs
    all_data_df = pd.concat([x.data_df for x in gctoos], axis=1)
    logger.info("all_data_df.shape: {}".format(all_data_df.shape))

    # Make sure df shapes are correct
    assert all_data_df.shape[0] == all_row_metadata_df.shape[0], "Number of rows is incorrect."
    assert all_data_df.shape[1] == all_col_metadata_df.shape[0], "Number of columns is incorrect."
    
    # Reset sample id index if requested (necessary if sample ids are not unique between plates)
    if reset_sample_ids:
        # See how many samples are repeated before resetting
        logger.info("num samples: {}".format(all_col_metadata_df.shape[0]))
        logger.info("num unique samples before reset: {}".format(len(np.unique(all_col_metadata_df.index.values))))

        # First, make sure sample ids agree between data df and col_meta_df
        assert(np.array_equal(all_col_metadata_df.index.values, all_data_df.columns.values))

        # Change index name so that the column that it becomes will be
        # appropriately named
        all_col_metadata_df.index.name = "old_cid"

        # Reset index
        all_col_metadata_df = all_col_metadata_df.reset_index()

        # Change the index name back to cid
        all_col_metadata_df.index.name = "cid"

        # Replace sample ids in data_df with the new ones from col_meta_df
        all_data_df.columns = pd.Index(all_col_metadata_df.index.values)

        # Confirm that num unique samples equals the number of samples
        logger.info("num unique samples after reset: {}".format(len(np.unique(all_col_metadata_df.index.values))))

        # Show the col_meta_df for debugging purposes
        logger.debug("all_col_metadata_df:\n{}".format(all_col_metadata_df))

    logger.info("build GCToo of all...")
    g = GCToo.GCToo(version="1.3", row_metadata_df=all_row_metadata_df, col_metadata_df=all_col_metadata_df, data_df=all_data_df)

    return g

def concat_row_meta(row_meta_dfs, fields_to_remove):
    """Concatenate the row metadata dfs together and sort the index.

    Args:
        row_meta_dfs (list of pandas dfs)
        fields_to_remove (list of strings): metadata fields to drop from all dfs
            before concatenating

    Returns:
        all_row_meta_df (pandas df)
    """
    # Remove any metadata columns that will prevent probes from being identical between plates
    if fields_to_remove is not None:
        for df in row_meta_dfs:
            df.drop(fields_to_remove, axis=1, inplace=True)

    # Concat all row_meta_df and then remove duplicate rows (slow...but it works)
    all_row_meta_df_dups = pd.concat(row_meta_dfs, axis=0)
    logger.debug("all_row_meta_df_dups.shape: {}".format(all_row_meta_df_dups.shape))
    all_row_meta_df = all_row_meta_df_dups.drop_duplicates()
    logger.debug("all_row_meta_df.shape: {}".format(all_row_meta_df.shape))

    # Verify that each rid now appears only once in the index
    assert all_row_meta_df.index.is_unique, (
        ("At least 1 rid has different metadata in different files, so it shows " +
         "up at least twice in all_row_meta_df. Correct this before proceeding. " +
         "all_row_meta_df:\n{}").format(all_row_meta_df))

    # Finally, re-sort the index
    all_row_metadata_df_sorted = all_row_meta_df.sort_index()

    return all_row_metadata_df_sorted

# def concat_col_meta(gctoos):
#     pass
#
# def concat_data(gctoos):
#     pass

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)
