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
import numpy as np
import os
import sys

import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.utils.psp_utils as psp_utils
import broadinstitute_psp.tear.continuous_renormalization as c
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.write_gct as wg

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

DEFAULT_TEAR_SUFFIX = ".tear.processed.gct"
ZSCORE_PROV_CODE_ENTRY = "ZSC"
PROV_CODE_FIELD = "provenance"
CONSTANT_FOR_MAD = 0.6745
MINIMUM_MAD = 1e-6


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--in_gct_path", "-i", required=True,
                        help="filepath to input gct")

    # Optional args
    parser.add_argument("--out_name", "-o", default=None,
                        help="name of output file (default is <INPUT_GCT>.tear.processed.gct")
    parser.add_argument("--divide_by_mad", "-dm", action="store_true", default=False,
                    help=("whether to divide by median absolute deviation " +
                          "in addition to subtracting the probe median"))
    parser.add_argument("--ignore_subset_norm", "-ig", action="store_true", default=False,
                        help="whether to ignore subset-specific normalization")
    parser.add_argument("-psp_config_path", type=str,
                        default="~/psp_production.cfg",
                        help="filepath to PSP config file")
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                        help="increase the number of messages reported")

    return parser


def main(args):
    # Read gct and config file
    (in_gct, config_io, config_metadata, _) = (
        psp_utils.read_gct_and_config_file(args.in_gct_path, args.psp_config_path))

    # Extract provenance code
    prov_code = psp_utils.extract_prov_code(
        in_gct.col_metadata_df, config_metadata["prov_code_field"],
        config_metadata["prov_code_delimiter"])

    ### MEDIAN NORMALIZE
    (out_gct, prov_code) = median_normalize(
        in_gct, args.divide_by_mad, args.ignore_subset_norm,
        config_metadata, prov_code)

    # Configure output name
    out_gct_name = configure_out_name(args.in_gct_path, args.out_name)

    # Reinsert provenance code
    out_gct.col_metadata_df = insert_prov_code(
        out_gct.col_metadata_df, prov_code,
        config_metadata["prov_code_delimiter"],
        config_metadata["prov_code_field"])

    # Apply continuous renormalization to GCT
    continuous_renormalization_args = c.build_parser(["--in_gct", out_gct])
    out_gct = c.continuous_renormalization(continuous_renormalization_args)

    # Write output gct
    write_output_gct(out_gct, out_gct_name, config_io["data_null"], config_io["filler_null"])
    return out_gct


# tested #
def median_normalize(gct, divide_by_mad, ignore_subset_norm, config_metadata, prov_code):
    """Subset normalize if the metadata shows that subsets exist for either
    rows or columns AND ignore_subset_norm is False. Otherwise, use the
    whole row for median normalization.

    Args:
        gct (GCToo object)
        divide_by_mad (bool): whether to divide by the median absolute deviation
            in addition to subtracting the median
        ignore_subset_norm (bool): false indicates that subset normalization should be performed
        config_metadata (dict): dictionary from config file with these fields:
            - row_subset_field
            - col_subset_field
            - subset_normalize_prov_code_entry
            - subset_zscore_prov_code_entry
            - row_normalize_prov_code_entry
            - zscore_prov_code_entry
        prov_code (list of strings)

    Returns:
        out_gct (GCToo object)
        updated_prov_code (list of strings)

    """
    if ignore_subset_norm:

        # Subtract median of whole row from each entry in the row
        out_df = row_median_normalize(gct.data_df, divide_by_mad)

        # Provenance code entry depends on whether z-scoring is happening
        if divide_by_mad:
            prov_code_entry = config_metadata["zscore_prov_code_entry"]
        else:
            prov_code_entry = config_metadata["row_normalize_prov_code_entry"]
        updated_prov_code = prov_code + [prov_code_entry]

    else:

        # Check if subsets_exist
        subsets_exist = check_for_subsets(gct.row_metadata_df, gct.col_metadata_df,
                                          config_metadata["row_subset_field"],
                                          config_metadata["col_subset_field"])
        if subsets_exist:
            logger.info("Subsets were found, so subset normalization will be performed.")

            out_df = subset_normalize(gct, divide_by_mad,
                                      config_metadata["row_subset_field"],
                                      config_metadata["col_subset_field"])

            # Provenance code entry depends on whether z-scoring is happening
            if divide_by_mad:
                prov_code_entry = config_metadata["subset_zscore_prov_code_entry"]
            else:
                prov_code_entry = config_metadata["subset_normalize_prov_code_entry"]
            updated_prov_code = prov_code + [prov_code_entry]

        else:
            logger.info("No subsets were found, so regular normalization will be performed.")

            # Subtract median of whole row from each entry in the row
            out_df = row_median_normalize(gct.data_df, divide_by_mad)

            # Provenance code entry depends on whether z-scoring is happening
            if divide_by_mad:
                prov_code_entry = config_metadata["zscore_prov_code_entry"]
            else:
                prov_code_entry = config_metadata["row_normalize_prov_code_entry"]
            updated_prov_code = prov_code + [prov_code_entry]


    out_gct = GCToo.GCToo(data_df=out_df, row_metadata_df=gct.row_metadata_df,
                          col_metadata_df=gct.col_metadata_df)

    return out_gct, updated_prov_code

# tested #
def check_for_subsets(row_meta_df, col_meta_df, row_subset_field, col_subset_field):
    """Read metadata to see if subset normalization could be performed. That is,
     if the row or column metadata has more than one unique value in
     the subset fields.

    Args:
        row_meta_df (pandas df)
        col_meta_df (pandas df)
        row_subset_field (string): row metadata field indicating the row subset group
        col_subset_field (string): col metadata field indicating the col subset group

    Returns:
        subsets_exist (bool): true if there

    """
    # Verify that subset metadata fields are in the metadata dfs
    assert row_subset_field in row_meta_df.columns, (
        "Row_meta_df does not have the row_subset_field: {}.".format(row_subset_field))
    assert col_subset_field in col_meta_df.columns, (
        "Col_meta_df does not have the col_subset_field: {}.".format(col_subset_field))

    num_unique_row_subsets = len(np.unique(row_meta_df.loc[:, row_subset_field].values));
    num_unique_col_subsets = len(np.unique(col_meta_df.loc[:, col_subset_field].values));

    # If either metadata field has more than one unique value, return true
    if num_unique_row_subsets > 1 or num_unique_col_subsets > 1:
        subsets_exist = True
    else:
        subsets_exist = False

    return subsets_exist


# tested #
def subset_normalize(gct, divide_by_mad, row_subset_field, col_subset_field):
    """Row-median normalize subsets of data.

    The function make_norm_ndarray extracts relevant metadata to figure out
    which subsets of the data should be row median normalized together, and
    iterate_over_norm_ndarray_and_normalize actually applies the normalization.

    Args:
        gct (GCToo object)
        divide_by_mad (bool): whether to divide by the median absolute deviation
            in addition to subtracting the median
        row_subset_field (string): row metadata field indicating the row subset group
        col_subset_field (string): col metadata field indicating the col subset group

    Returns:
        out_df (pandas df): with normalized data
    """
    # Create normalization ndarray
    norm_ndarray = make_norm_ndarray(gct.row_metadata_df, gct.col_metadata_df, row_subset_field, col_subset_field)

    # Iterate over norm ndarray and actually perform normalization
    out_df = iterate_over_norm_ndarray_and_normalize(gct.data_df, norm_ndarray, divide_by_mad)

    return out_df

# tested #
def make_norm_ndarray(row_metadata_df, col_metadata_df, row_subset_field, col_subset_field):
    """Creates a normalization ndarray.

    The normalization ndarray is used to figure out how to normalize each
    section of the data. Extracts probe normalization groups from the metadata
    field row_subset_field and from the column metadata field col_subset_field.

    Args:
        row_metadata_df (pandas df)
        col_metadata_df (pandas df)
        row_subset_field (string): row metadata field indicating the row subset group
        col_subset_field (string): col metadata field indicating the col subset group

    Returns:
        norm_ndarray (numpy array): size = (num_probes, num_samples)
    """
    # Verfiy that metadata fields of interest exist
    assert row_subset_field in row_metadata_df.columns
    assert col_subset_field in col_metadata_df.columns

    # Get sample group vectors
    sample_grps_strs = col_metadata_df[col_subset_field].values

    # Vectors could be strings, or could have been converted to integers if
    # the vector had length 1; so convert to str just in case
    sample_grps_strs = sample_grps_strs.astype(str)

    # Convert sample group vectors from strs to lists of strings
    sample_grps_lists = [sample_str.split(",") for sample_str in sample_grps_strs]

    # Verify that all lists have same length
    length_of_first_list = len(sample_grps_lists[0])
    assert all([len(sample_list) == length_of_first_list for sample_list in sample_grps_lists])

    # Convert from lists of strings to ndarray of ints
    sample_grp_ndarray = np.array(sample_grps_lists, dtype="float").astype("int")

    # Get probe groups and unique probe groups; convert to ints
    probe_grps = row_metadata_df[row_subset_field].values.astype("int")
    unique_probe_grps = np.unique(probe_grps)

    # Length of sample group vector must equal number of unique probe groups
    assert length_of_first_list == len(unique_probe_grps), (
        "Length of sample group vector must equal number of unique " +
        "probe groups. len(sample_grps_lists[0]): {} \n " +
        "len(unique_probe_grps): {}".format(len(sample_grps_lists[0]), len(unique_probe_grps)))

    # Initialize norm_ndarray
    norm_ndarray = np.zeros((row_metadata_df.shape[0], col_metadata_df.shape[0]), dtype="int")

    # Each col of sample_grp_ndarray corresponds to a UNIQUE probe_grp;
    # create norm_ndarray by expanding sample_grp_ndarray so that each column
    # corresponds to a probe
    for unique_probe_grp_idx in range(len(unique_probe_grps)):
        # Row to insert is the column from sample_grp_ndarray
        row_to_insert = sample_grp_ndarray[:, unique_probe_grp_idx]

        # Select probe indices for which this row should be inserted
        probe_bool_vec = (probe_grps == unique_probe_grps[unique_probe_grp_idx])

        # Insert into norm_ndarray
        norm_ndarray[probe_bool_vec, :] = row_to_insert

    # Check that there are no zeros in norm_ndarray anymore
    assert ~np.any(norm_ndarray == 0), (
        "There should not be any zeros in norm_ndarray anymore.")

    return norm_ndarray

# tested #
def iterate_over_norm_ndarray_and_normalize(data_df, norm_ndarray, divide_by_mad):
    """Iterate over norm_ndarray and row-median normalize the subsets.

    Each row of norm_ndarray indicates the subsets for normalization to
    select from data_df.

    Args:
        data_df (pandas df): size = (num_probes, num_samples)
        norm_ndarray (numpy ndarray): size = (num_probes, num_samples)
        divide_by_mad (bool): whether to divide by the median absolute deviation
            in addition to subtracting the median
    Returns:
        normalized_data (pandas df): values will be modified
    """
    normalized_data = data_df.copy()

    for row_idx, norm_ndarray_row in enumerate(norm_ndarray):
        # Get number of subsets for this probe
        sample_subsets = np.unique(norm_ndarray_row)

        for sample_subset in sample_subsets:
            # Select data to normalize
            data_to_normalize = data_df.iloc[row_idx].loc[norm_ndarray_row == sample_subset].values

            # Subtract median and divide by MAD
            if divide_by_mad:
                mad = np.nanmedian(np.absolute(data_to_normalize - np.nanmedian(data_to_normalize)))

                # Set some minimum value for the MAD so as not to divide by 0
                mad = max(mad, MINIMUM_MAD)

                data_to_insert = (data_to_normalize - np.nanmedian(data_to_normalize)) / (CONSTANT_FOR_MAD * mad)

            # Simply subtract median
            else:
                data_to_insert = data_to_normalize - np.nanmedian(data_to_normalize)

            # Insert data into output ndarray
            normalized_data.iloc[row_idx].loc[norm_ndarray_row == sample_subset] = data_to_insert

    return normalized_data

# tested #
def row_median_normalize(data_df, divide_by_mad):
    """Subtract median of the row from each entry of the row. If divide_by_mad=True,
    also divide by the median absolute deviation.

    Args:
        data_df (pandas df)
        divide_by_mad (bool): whether to divide by the median absolute deviation
            in addition to subtracting the median
    Returns:
        out_df (pandas df): with normalized values
    """
    if divide_by_mad:

        # Calculate MAD
        mad = np.absolute(data_df.subtract(data_df.median(axis=1), axis=0)).median(axis=1)

        # Set some minimum value for the MAD so as not to divide by 0
        mad = mad.apply(max, args=(MINIMUM_MAD,))

        # Subtract median and divide by MAD
        out_df = (data_df.subtract(data_df.median(axis=1), axis=0)).divide(CONSTANT_FOR_MAD * mad, axis=0)

    # Simply subtract median
    else:
        out_df = data_df.subtract(data_df.median(axis=1), axis=0)

    return out_df


def insert_prov_code(col_metadata_df, prov_code, prov_code_delimiter, prov_code_field):
    """ Update provenance code in col_metadata_df.

    Args:
        col_metadata_df (pandas df)
        prov_code (list of strings)
        prov_code_field (string): name of col metadata field containing the provenance code
        prov_code_delimiter (string): what string to use as delimiter in prov_code

    Returns:
        gct (GCToo object): updated metadata

    """
    # Convert provenance code to delimiter separated string
    prov_code_str = prov_code_delimiter.join(prov_code)

    # Update the provenance code in col_metadata_df
    col_metadata_df.loc[:, prov_code_field] = prov_code_str

    return col_metadata_df


def configure_out_name(in_gct_path, out_name_from_args):
    """If out_name_from_args is None, append DEFAULT_TEAR_SUFFIX to the input
    gct name.

    Args:
        in_gct_path (file path)
        out_name_from_args (string)

    Returns:
        out_gct_name (file path)

    """
    input_basename = os.path.basename(in_gct_path)

    if out_name_from_args is None:
        out_gct_name = input_basename + DEFAULT_TEAR_SUFFIX
    else:
        out_gct_name = out_name_from_args
        assert os.path.splitext(out_gct_name)[1] == ".gct", (
            "The output gct name must end with .gct; out_gct_name: {}".format(
                out_gct_name))

    return out_gct_name


def write_output_gct(out_gct, out_gct_name, data_null, filler_null):

    wg.write(out_gct, out_gct_name, data_null=data_null, filler_null=filler_null, data_float_format=None)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args: {}".format(args))

    main(args)