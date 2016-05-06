import logging
import utils.setup_logger as setup_logger
import sys
import argparse
import in_out.parse_gctoo as parse_gctoo
import pandas as pd

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

"""
Converts a QC gct file into a .pw (plate-well) file that can be used
easily with Plato.

Divides each value by the maximum value in the row, and then computes the
median, mean, MAD, and SD of the values in each column.

Input is a gct file. Output is a pw file.
"""

# TODO(lev): update documentation to reflect correct computations

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

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
    gct = parse_gctoo.parse(args.gct_file_path)

    # Get plate and well names
    (plate_names, well_names) = extract_plate_and_well_names(
        gct.col_metadata_df, args.plate_field, args.well_field)

    # Divide by the maximum value for the row
    max_row_values = gct.data_df.max(axis='columns')
    divided_data_df = gct.data_df.div(max_row_values, axis="rows")

    # TODO(lev): what do I do with negative values?

    # Calculate metrics for each sample
    medium_over_heavy_medians = divided_data_df.median(axis=0).values
    medium_over_heavy_means = divided_data_df.mean(axis=0).values
    medium_over_heavy_mads = divided_data_df.mad(axis=0).values
    medium_over_heavy_sds = divided_data_df.std(axis=0).values

    # Assemble plate_names, well_names, and metrics into a dataframe
    out_df = assemble_output_df(
        plate_names, well_names,
        medium_over_heavy_median=medium_over_heavy_medians,
        medium_over_heavy_mad=medium_over_heavy_mads)

    # Write to pw file
    out_df.to_csv(args.out_pw_file_path, sep="\t", na_rep="NaN", index=False)
    logger.info("PW file written to {}".format(args.out_pw_file_path))

def extract_plate_and_well_names(col_meta, plate_field, well_field):
    """

    Args:
        col_meta (pandas df)
        plate_field (string):
        well_field (string):

    Returns:
        plate_names (numpy array of strings)
        well_names (numpy array of strings)

    """

    # Extract plate metadata
    plate_names = col_meta[plate_field].values

    # Verify that all samples have the same plate name
    plate_names_same = True
    for plate in plate_names:
        plate_names_same = (plate_names_same and plate == plate_names[0])
        assert plate_names_same

    # Extract well metadata
    well_names = col_meta[well_field].values

    return plate_names, well_names


def assemble_output_df(plate_names, well_names, **kwargs):
    """Assemble output df for saving.

    plate_name will be the first column, well_name will be the second column,
    and the remaining columns will be alphabetically ordered by kwargs.keys().

    Args:
        plate_names (numpy array of strings)
        well_names (numpy array of strings)
        kwargs (dict of keyword pairs): values must be iterables with length = num wells

    Returns:
        out_df (pandas df)
    """
    PLATE_FIELD = "plate_name"
    WELL_FIELD = "well_name"

    # Make plate and well names into a dict
    plate_and_well_dict = {PLATE_FIELD: plate_names, WELL_FIELD: well_names}

    assert len(plate_names) == len(well_names)
    num_wells = len(well_names)

    # Assert that length of each optional value is equal to number of wells
    for kwarg in kwargs.items():
        assert len(kwarg[1]) == num_wells, (
            "The entry {} has length {}, which is not equal to num_wells: {}.".format(
                kwarg[0], len(kwarg[1]), num_wells))

    # Append the kwargs dict to plate_and_well_dict
    df_dict = plate_and_well_dict.copy()
    df_dict.update(kwargs)

    # Convert dict to df and rearrange columns appropriately
    temp_df = pd.DataFrame.from_dict(df_dict)
    cols = [PLATE_FIELD, WELL_FIELD] + sorted(kwargs.keys())
    out_df = temp_df[cols]

    return out_df


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup()
    main(args)