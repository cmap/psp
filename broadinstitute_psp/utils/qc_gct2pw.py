"""
QC_GCT2PW.PY: Converts a QC gct file into a .pw (plate-well) file that can
be used easily with Plato.

Divides each value by the maximum value of its row, and then computes the
median, mean, MAD, and SD for each column.

Input is a gct file. Output is a pw file.
"""

import logging
import sys
import argparse
import numpy as np
import pandas as pd

import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.utils.psp_utils as utils
import cmapPy.pandasGEXpress.parse as parse

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

PLATE_FIELD = "plate_name"
WELL_FIELD = "well_name"
PROV_CODE_FIELD = "provenance_code"
PROV_CODE_DELIMITER = "+"
LOG_TRANSFORM_PROV_CODE_ENTRY = "L2X"


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("gct_file_path", type=str, help="filepath to gct file")
    parser.add_argument("out_pw_file_path", type=str,
                        help="filepath to output pw file")
    
    # Optional args
    parser.add_argument("-plate_field", type=str, default="det_plate",
                        help="metadata field name specifying the plate")
    parser.add_argument("-well_field", type=str, default="det_well",
                        help="metadata field name specifying the well")
    return parser


def main(args):
    # Import gct
    gct = parse.parse(args.gct_file_path)

    # Get plate and well names
    (plate_names, well_names) = extract_plate_and_well_names(
        gct.col_metadata_df, args.plate_field, args.well_field)

    # Extract provenance code
    prov_code = utils.extract_prov_code(
        gct.col_metadata_df, PROV_CODE_FIELD, PROV_CODE_DELIMITER)

    # If data has been log-transformed, undo it
    unlogged_df = undo_log_transform_if_needed(gct.data_df, prov_code)

    # Divide by the maximum value for the row
    max_row_values = unlogged_df.max(axis='columns')
    divided_df = unlogged_df.div(max_row_values, axis="rows")

    # Calculate metrics for each sample
    medium_over_heavy_medians = divided_df.median(axis=0).values
    medium_over_heavy_means = divided_df.mean(axis=0).values
    medium_over_heavy_mads = divided_df.mad(axis=0).values
    medium_over_heavy_sds = divided_df.std(axis=0).values

    # Assemble plate_names, well_names, and metrics into a dataframe
    out_df = assemble_output_df(
        plate_names, well_names,
        {"medium_over_heavy_median": medium_over_heavy_medians,
        "medium_over_heavy_mad": medium_over_heavy_mads})

    # Write to pw file
    out_df.to_csv(args.out_pw_file_path, sep="\t", na_rep="NaN", index=False)
    logger.info("PW file written to {}".format(args.out_pw_file_path))

def extract_plate_and_well_names(col_meta, plate_field=PLATE_FIELD, well_field=WELL_FIELD):
    """

    Args:
        col_meta (pandas df)
        plate_field (string): metadata field for name of plate
        well_field (string): metadata field for name of well

    Returns:
        plate_names (numpy array of strings)
        well_names (numpy array of strings)

    """

    # Extract plate metadata
    plate_names = col_meta[plate_field].values

    # Extract well metadata
    well_names = col_meta[well_field].values

    return plate_names, well_names


def assemble_output_df(plate_names, well_names, metadata_dict):
    """Assemble output df for saving.

    plate_name will be the first column, well_name will be the second column,
    and the remaining columns will be alphabetically ordered by
    extra_metadata_dict.keys().

    Args:
        plate_names (numpy array of strings)
        well_names (numpy array of strings)
        metadata_dict (dictionary): keys will become the names of metadata
            fields, and values must be iterables with length = num wells

    Returns:
        out_df (pandas df)
    """

    # Make plate and well names into a dict
    plate_and_well_dict = {PLATE_FIELD: plate_names, WELL_FIELD: well_names}

    assert len(plate_names) == len(well_names)
    num_wells = len(well_names)

    # Assert that length of each metadata value is equal to number of wells
    for meta_key, meta_value in metadata_dict.iteritems():
        assert len(meta_value) == num_wells, (
            "The entry {} has length {}, which is not equal to num_wells: {}.".format(
                meta_key, len(meta_value), num_wells))

    # Append the metadata_dict to plate_and_well_dict
    df_dict = plate_and_well_dict.copy()
    df_dict.update(metadata_dict)

    # Convert dict to df and rearrange columns appropriately
    temp_df = pd.DataFrame.from_dict(df_dict)
    cols = [PLATE_FIELD, WELL_FIELD] + sorted(metadata_dict.keys())
    out_df = temp_df[cols]

    return out_df

def undo_log_transform_if_needed(data_df, prov_code):
    """Undo log transformation if L2X is in prov_code."""
    if LOG_TRANSFORM_PROV_CODE_ENTRY in prov_code:
        out_df = np.exp2(data_df)
    else:
        out_df = data_df

    return out_df


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup()
    main(args)
