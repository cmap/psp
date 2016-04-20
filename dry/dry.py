import logging
import utils.setup_logger as setup_logger
import ConfigParser
import argparse
import sys
import os
import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import in_out.GCToo as GCToo
import in_out.parse_gctoo as parse_gctoo
import in_out.write_gctoo as write_gctoo

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

"""Performs filtering and normalization of P100 and GCP data.

Converts level 2 to level 3 data. Required input is a filepath to a gct file,
a filepath to an output file, and what to name the output file.
Output is writing the gct file.

N.B. This script requires a configuration file. You can specify the location
of this config file with the optional argument -psp_config_path. Otherwise,
it will look for psp.cfg in your home directory. See example_psp.cfg for
an example config file.

Example usage:
    python dry/dry.py input.gct /path/to/output/dir output.gct
    -sample_nan_thresh 0.9 -subset_normalize_bool -optim -force_assay PR1

"""
# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

LOG_TRANSFORM_PROV_CODE_ENTRY = "L2X"
GCP_HISTONE_PROV_CODE_ENTRY = "H3N"
SAMPLE_FILTER_PROV_CODE_ENTRY = "SF"
MANUAL_PROBE_FILTER_PROV_CODE_ENTRY = "MPR" # not used
PROBE_FILTER_PROV_CODE_ENTRY = "PF"
OPTIMIZATION_PROV_CODE_ENTRY = "LLB" # rename to OPT?
OUTLIER_SAMPLE_FILTER_PROV_CODE_ENTRY = "OSF" # add sd cutoff?
SUBSET_NORMALIZE_PROV_CODE_ENTRY = "GMN" # rename to SSN (SubSetNormalize)?
ROW_NORMALIZE_PROV_CODE_ENTRY = "RMN"

def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("gct_path", type=str, help="filepath to gct file")
    parser.add_argument("out_path", type=str, help="path to save directory")
    parser.add_argument("out_name", type=str, help="name of output gct")

    # Optional args
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                        help="increase the number of messages reported")
    parser.add_argument("-psp_config_path", type=str,
                        default="~/psp.cfg",
                        help="filepath to PSP config file")
    parser.add_argument("-force_assay", choices=["GR1", "PR1"], default=None, # TODO(lev): should I allow other options?
                        help=("directly specify assay type here " +
                              "(overrides first entry of provenance code)"))
    parser.add_argument("-optim", action="store_true", default=False,
                        help="whether to perform load balancing optimization")
    parser.add_argument("-subset_normalize_bool", action="store_true", default=False,
                        help="whether to perform subset-specific normalization")
    # TODO(lev): do we want to keep this?

    # Parameters
    parser.add_argument("-sample_nan_thresh", type=float, default=0.8,
                        help=("if < sample_nan_thresh of a sample's data " +
                              "is non-nan, that sample will be filtered out"))
    parser.add_argument("-probe_nan_thresh", type=float, default=0.9,
                        help=("if < probe_nan_thresh of a probe's data " +
                              "is non-nan, that probe will be filtered out"))
    parser.add_argument("-probe_sd_cutoff", type=float, default=3,
                        help=("maximum SD for a probe before being filtered out"))
    parser.add_argument("-dist_sd_cutoff", type=float, default=3,
                        help=("maximum SD for a sample's distance metric " +
                              "before being filtered out"))
    parser.add_argument("-optim_bounds", type=tuple, default=(-7,7),
                        help="bounds over which to perform optimization")

    return parser


def main(args):
    """THE MAIN METHOD. Filter and normalize data and save result as a gct file.

    Args:
        args (argparse.Namespace object): fields as defined in build_parser()

    Returns:
        out_gct (GCToo object): output gct object
    """
    ### READ GCT AND CONFIG FILE
    (gct, assay_type, prov_code, p100_assay_types,
     gcp_assay_types, config_io, config_metadata) = (
        read_gct_and_config_file(
            args.gct_path, args.psp_config_path, args.force_assay))

    ### LOG TRANSFORM
    (gct.data_df, prov_code) = log_transform_if_needed(gct.data_df, prov_code)

    ### HISTONE NORMALIZE (if GCP)
    (gct, prov_code) = gcp_histone_normalize_if_needed(
        gct, assay_type, gcp_assay_types,
        config_metadata["gcp_normalization_peptide_id"], prov_code)

    ### INITIAL FILTERING
    (gct, prov_code) = initial_filtering(
        gct, args.sample_nan_thresh, args.probe_nan_thresh,
        args.probe_sd_cutoff, config_metadata["manual_rejection_field"], prov_code)

    ### APPLY OFFSETS IF NEEDED (if P100)
    (gct, dists, offsets, success_bools, prov_code) = (
        p100_calculate_dists_and_apply_offsets_if_needed(
            gct, assay_type, p100_assay_types, args.optim,
            args.optim_bounds, prov_code))

    ### FILTER SAMPLES BY DISTANCE (if P100)
    (gct, out_offsets, prov_code) = p100_filter_samples_by_dist(
        gct, assay_type, p100_assay_types, offsets, dists,
        success_bools, args.dist_sd_cutoff, prov_code)

    ### MEDIAN NORMALIZE
    (gct, prov_code) = median_normalize(
        gct, args.subset_normalize_bool, config_metadata["row_subset_field"],
        config_metadata["col_subset_field"], prov_code)

    ### INSERT OFFSETS AND UPDATE PROVENANCE CODE
    out_gct = insert_offsets_and_prov_code(
        gct, out_offsets, config_metadata["offsets_field"], prov_code,
        config_metadata["prov_code_field"], config_metadata["prov_code_delimiter"])

    ### WRITE OUTPUT GCT
    write_output_gct(out_gct, args.out_path, args.out_name,
                     config_io["data_null"], config_io["filler_null"])

    return out_gct


# tested #
def read_gct_and_config_file(gct_path, config_path, forced_assay_type):
    """Read gct and config file.

    The config file has two sections: io and metadata. These are returned as
    dictionaries. The field "nan_values" in "io" indicates what values to
    consider NaN when reading in the gct file. The fields "gcp_assays" and
    "p100_assays" are returned in order to be used by later functions that need
    to know if the assay is P100 or GCP.

    Provenance code is extracted from the col metadata. It must be non-empty
    and the same for all samples. If forced_assay_type is not None,
    assay_type is set to forced_assay_type.

    Args:
        gct_path (string): filepath to gct file
        config_path (srting): filepath to config file
        forced_assay_type (string, or None)
        prov_code_field (string): name of the col metadata field for the provenance code

    Returns:
        gct (GCToo object)
        assay_type (string)
        prov_code (list of strings)
        p100_assay_types (list of strings)
        gcp_assay_types (list of strings)
        config_io (dictionary)
        config_metadata (dictionary)
    """

    # Read from config file
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(os.path.expanduser(config_path))
    config_io = dict(config_parser.items("io"))
    config_metadata = dict(config_parser.items("metadata"))

    # Extract what values to consider NaN
    # N.B. eval used to convert string from config file to list
    psp_nan_values = eval(config_io["nan_values"])

    # Parse the gct file and return GCToo object
    gct = parse_gctoo.parse(gct_path, nan_values=psp_nan_values)

    # Extract the plate's provenance code
    prov_code = extract_prov_code(gct.col_metadata_df,
                                  config_metadata["prov_code_field"],
                                  config_metadata["prov_code_delimiter"])

    # If forced_assay_type is not None, set assay_type to forced_assay_type.
    # Otherwise, the first entry of the provenance code is the assay_type.
    if forced_assay_type is not None:
        assay_type = forced_assay_type
    else:
        assay_type = prov_code[0]

    # Make sure assay_type is one of the allowed values
    p100_assay_types = eval(config_metadata["p100_assays"])
    gcp_assay_types = eval(config_metadata["gcp_assays"])
    allowed_assay_types = p100_assay_types + gcp_assay_types
    check_assay_type(assay_type, allowed_assay_types)

    return gct, assay_type, prov_code, p100_assay_types, gcp_assay_types, config_io, config_metadata


# tested #
def extract_prov_code(col_metadata_df, prov_code_field, prov_code_delimiter):
    """Extract the provenance code from the column metadata.
    It must be non-empty and the same for all samples.

    Args:
        col_metadata_df (pandas df): contains provenance code metadata
        prov_code_field (string): name of the col metadata field for the provenance code
        prov_code_delimiter (string): string delimiter in provenance code

    Returns:
        prov_code (list of strings)

    """
    # Create pandas series of all provenance codes
    prov_code_series = col_metadata_df.loc[:, prov_code_field]

    # Split each provenance code string along the delimiter
    prov_code_list_series = prov_code_series.apply(lambda x: x.split(prov_code_delimiter))

    # Verify that all provenance codes are the same
    # (i.e. verify that the number of elements equal to the first element is
    # equal to the number of all elements)
    all_same = True
    for prov_code_list in prov_code_list_series:
        all_same = (all_same and prov_code_list == prov_code_list_series[0])

    if all_same:
        prov_code = prov_code_list_series.iloc[0]
        assert prov_code != [""], "Provenance code is empty!"
        return prov_code
    else:
        err_msg = ("All columns should have the same provenance code, " +
                   "but actually np.unique(prov_code_list_series.values) = {}")
        logger.error(err_msg.format(np.unique(prov_code_list_series.values)))
        raise(Exception(err_msg.format(np.unique(prov_code_list_series.values))))


# tested #
def check_assay_type(assay_type, allowed_assay_types):
    """Verify that assay is one of the allowed types.

    Args:
        assay_type (string)
        allowed_assay_types (list of strings)

    Returns:
        assay_ok (bool)
    """
    assay_ok = (assay_type in allowed_assay_types)
    if not assay_ok:
        err_msg = ("The assay type is not in the allowed assay types. " +
                   "assay_type: {}, allowed_assay_types: {}")
        logger.error(err_msg.format(assay_type, allowed_assay_types))
        raise(Exception(err_msg.format(assay_type, allowed_assay_types)))
    return assay_ok


# tested #
def log_transform_if_needed(data_df, prov_code):
    """Perform log2 transformation if it hasn't already been done.

    Args:
        data_df (pandas df)
        prov_code (list of strings)

    Returns:
        out_df (pandas df)
        updated_prov_code (list of strings): updated
    """
    # Check if log2 transformation has already occurred
    if LOG_TRANSFORM_PROV_CODE_ENTRY in prov_code:
        logger.info("{} has already occurred.".format(LOG_TRANSFORM_PROV_CODE_ENTRY))
        out_df = data_df
        updated_prov_code = prov_code
    else:
        out_df = log_transform(data_df, log_base=2)
        prov_code_entry = LOG_TRANSFORM_PROV_CODE_ENTRY
        updated_prov_code = prov_code + [prov_code_entry]

    return out_df, updated_prov_code


# tested #
def log_transform(data_df, log_base):
    """Take the log of each value in a pandas dataframe.

    Args:
        data_df (pandas df)
        log_base (int): the base of the logarithm
    Returns:
        out_df (pandas df)
    """
    # Replace 0 with np.nan
    data_df_no_zeros = data_df.replace(0, np.nan)

    # Numpy operations work fine on dataframes
    out_df = np.log(data_df_no_zeros) / np.log(log_base)

    return out_df

# tested #
def gcp_histone_normalize_if_needed(gct, assay_type, gcp_assays, gcp_normalization_peptide_id, prov_code):
    """If GCP and gcp_normalization_peptide_id is not None, perform
     histone normalization. This subtracts the normalization probe row
     from all the other probe rows. It also removes the normalization probe
     from the data going forward.

    Args:
        gct (GCToo object)
        assay_type (string)
        gcp_assays (list of strings)
        gcp_normalization_peptide_id (string, or None)
        prov_code (list of strings)

    Returns:
        out_gct (GCToo object): updated data and metadata dfs;
            one row removed from each if normalization occurred
        updated_prov_code (list of strings)
    """
    if (assay_type in gcp_assays) and (gcp_normalization_peptide_id is not None):
        # Perform normalization
        out_df = gcp_histone_normalize(gct.data_df, gcp_normalization_peptide_id)

        # Update gct object and provenance code
        prov_code_entry = GCP_HISTONE_PROV_CODE_ENTRY
        (out_gct, updated_prov_code) = update_metadata_and_prov_code(
            out_df, gct.row_metadata_df, gct.col_metadata_df, prov_code_entry, prov_code)

    else:
        out_gct = gct
        updated_prov_code = prov_code

    return out_gct, updated_prov_code


# tested #
def gcp_histone_normalize(data_df, gcp_normalization_peptide_id):
    """Subtract values of gcp_normalization_peptide_id from all the other probes.
    Remove the row of data corresponding to the normalization histone.

    Args:
        data_df (pandas df)
        gcp_normalization_peptide_id (string)

    Returns:
        out_df (pandas df): one row removed
    """
    # Verify that normalization peptide is in the data
    assert gcp_normalization_peptide_id in data_df.index, (
        ("The normalization peptide is not in this dataset. " +
         "gcp_normalization_peptide_id: {}".format(gcp_normalization_peptide_id)))

    # Calculate normalization values
    norm_values = data_df.loc[gcp_normalization_peptide_id, :]

    # Drop the normalization peptide row
    out_df = data_df.drop(gcp_normalization_peptide_id)

    # Subtract the normalization values from all rows
    out_df = out_df - norm_values
    return out_df

# tested #
def initial_filtering(gct, sample_nan_thresh, probe_nan_thresh, probe_sd_cutoff, manual_rejection_field, prov_code):
    """Performs three types of filtering. filter_samples_by_nan removes
    samples that have too many nan values. manual_probe_rejection removes
    probes that were manually labeled for removal. filter_probes_by_nan_and_sd
    removes probes if they have too many nan values OR if the SD for the whole
    row is too high.

    Args:
        gct (GCToo object)
        sample_nan_thresh (float b/w 0 and 1)
        probe_nan_thresh (float b/w 0 and 1)
        probe_sd_cutoff (float)
        manual_rejection_field (string)
        prov_code (list of strings)

    Returns:
        gct (GCToo object): updated
        prov_code (list of strings): updated
    """

    # Filter samples by nan
    gct.data_df = filter_samples_by_nan(gct.data_df, sample_nan_thresh)
    thresh_digit = ("{:.1f}".format(sample_nan_thresh)).split(".")[1]
    prov_code_entry = "{}{}".format(SAMPLE_FILTER_PROV_CODE_ENTRY, thresh_digit)
    (gct, prov_code) = update_metadata_and_prov_code(
        gct.data_df, gct.row_metadata_df, gct.col_metadata_df, prov_code_entry, prov_code)

    # Filter manually rejected probes
    gct.data_df = manual_probe_rejection(gct.data_df, gct.row_metadata_df, manual_rejection_field)
    prov_code_entry = None # TODO(lev): do we want MPR as a provenance code entry?
    (gct, prov_code) = update_metadata_and_prov_code(
        gct.data_df, gct.row_metadata_df, gct.col_metadata_df, prov_code_entry, prov_code)

    # TODO(lev): should we separate filtering by nan and by sd?
    # Filter probes by nan and sd
    gct.data_df = filter_probes_by_nan_and_sd(gct.data_df, probe_nan_thresh, probe_sd_cutoff)
    thresh_digit = ("{:.1f}".format(probe_nan_thresh)).split(".")[1]
    prov_code_entry = "{}{}".format(PROBE_FILTER_PROV_CODE_ENTRY, thresh_digit)
    (gct, prov_code) = update_metadata_and_prov_code(
        gct.data_df, gct.row_metadata_df, gct.col_metadata_df, prov_code_entry, prov_code)

    return gct, prov_code


# tested #
def filter_samples_by_nan(data_df, sample_nan_thresh):
    """Remove samples (i.e. columns) with less than sample_nan_thresh non-NaN values.

    Args:
        data_df (pandas df)
        sample_nan_thresh (float b/w 0 and 1)
    Returns:
        out_df (pandas df): potentially smaller than input
    """
    # Number of NaNs per sample
    num_nans = data_df.isnull().sum()

    # Number of rows
    num_rows = data_df.shape[0]

    # Percent non-NaN per sample
    pct_non_nan_per_sample = 1 - num_nans / num_rows

    # Only return samples with fewer % of NaNs than the cutoff
    out_df = data_df.loc[:, pct_non_nan_per_sample > sample_nan_thresh]
    assert not out_df.empty, "All samples were filtered out. Try reducing the threshold."
    return out_df


# tested #
def manual_probe_rejection(data_df, row_metadata_df, manual_rejection_field):
    """Remove probes that were manually selected for rejection.

    Args:
        data_df (pandas df)
        row_metadata_df (pandas df)
        manual_rejection_field (string): row metadata field
    Returns:
        out_df (pandas df): potentially smaller than input
    """
    # Extract manual_rejection_field
    keep_probe_str = row_metadata_df.loc[:, manual_rejection_field]

    # Convert strings to booleans
    keep_probe_bool = (keep_probe_str == "TRUE")

    # Check that the list of good probes is not empty
    assert keep_probe_bool.any(), (
        "No probes were marked TRUE (i.e. suitable).\n" +
        "row_metadata_df.loc[:, '{}']: \n{}").format(manual_rejection_field, keep_probe_str)

    # Return the good probes
    out_df = data_df[keep_probe_bool.values]
    assert not out_df.empty, "All probes were filtered out!"
    return out_df

# tested #
def filter_probes_by_nan_and_sd(data_df, probe_nan_thresh, probe_sd_cutoff):
    """Remove probes (i.e. rows) with less than probe_nan_thresh non-NaN values.
    Also remove probes with standard deviation higher than probe_sd_cutoff.

    Args:
        data_df (pandas df)
        probe_nan_thresh (float b/w 0 and 1)
        probe_sd_cutoff (float)
    Returns:
        out_df (pandas df): potentially smaller than original df
    """
    # Input should be a pandas dataframe
    assert isinstance(data_df, pd.DataFrame), (
        "data_df must be a pandas dataframe, not type(data_df): {}").format(type(data_df))

    # Number of NaNs per probe
    num_nans = data_df.isnull().sum(axis=1)

    # Number of samples
    num_samples = data_df.shape[1]

    # Percent non-NaN per probe
    pct_non_nans_per_probe = 1 - num_nans / num_samples

    # Probe standard deviations
    probe_sds = data_df.std(axis=1)

    # Only return probes with fewer % of NaNs than the threshold
    # and lower sd than the cutoff
    probes_to_keep = ((pct_non_nans_per_probe > probe_nan_thresh) &
                      (probe_sds < probe_sd_cutoff))
    out_df = data_df.loc[probes_to_keep, :]
    assert not out_df.empty, (
        "All probes were filtered out. Try reducing the NaN threshold and/or SD cutoff.")
    return out_df

# tested #
def p100_calculate_dists_and_apply_offsets_if_needed(gct, assay_type,
                                                     p100_assays, optim_bool,
                                                     optim_bounds, prov_code):
    """If P100, calculate a distance metric for each sample in order to use it
    later for filtering. The distance metric is how far away each probe is from
    its median value.

    If optimization is also requested (optim_bool=True), an offset is calculated
    for each sample that seeks to minimize the distance metric for that sample.
    The distance is then recalculated for the data after offsets have been
    applied. The provenance code is only updated if optimization occurs.

    Args:
        gct (GCToo object)
        assay_type (string)
        p100_assays (list of strings): list of assay types considered P100
        optim_bool (bool): whether to perform optimization
        optim_bounds (tuple of floats): the bounds over which optimization should be performed
        prov_code (list of strings)

    Returns:
        gct (GCToo object): updated
        dists (numpy array of floats): distance metric for each sample
        offsets (numpy array of floats): offset that was applied to each sample
        success_bools (numpy array of bools): for each sample, indicates whether optimization was successful
        prov_code (list of strings): updated

    """
    # P100
    if assay_type in p100_assays:
        if optim_bool:
            # Perform optimization and return offsets, distances, and success_bools
            (gct.data_df, offsets, dists, success_bools) = (
                calculate_distances_and_optimize(gct.data_df, optim_bounds))
            prov_code_entry = OPTIMIZATION_PROV_CODE_ENTRY
            (gct, prov_code) = update_metadata_and_prov_code(
                gct.data_df, gct.row_metadata_df, gct.col_metadata_df, prov_code_entry, prov_code)
        else:
            # Simply calculate distance metric for each sample
            dists = calculate_distances_only(gct.data_df)

            # Set offsets and success_bools to None
            offsets = None
            success_bools = None

    # GCP
    # N.B. distances are not calculated because filtration by distance doesn't occur
    else:
        dists = None
        offsets = None
        success_bools = None

    return gct, dists, offsets, success_bools, prov_code


# tested #
def calculate_distances_and_optimize(data_df, optim_bounds):
    """For each sample, perform optimization.

    This means finding an offset for each sample that minimizes
    the distance metric. The distance metric is the sum of the
    distances from each probe measurement to its median.

    N.B. Only uses the non-NaN values from each sample.

    Args:
        data_df (pandas df)
        optim_bounds (tuple of floats): indicates range of values over which to perform optimization

    Returns:
        out_df (pandas dataframe): with offsets applied
        optimized_offsets (numpy array): constant to add to each sample to
            minimize the distance function; length = num_samples
        optim_distances (numpy array): distance with optimized offset added;
            length = num_samples
        success_bools: numpy array of bools indicating if optimization converged
            with length = num_samples
    """
    # Determine the median value for each probe
    probe_medians = data_df.median(axis=1)

    # Initialize optimization outputs
    num_samples = data_df.shape[1]
    optimized_offsets = np.zeros(num_samples, dtype=float)
    optimized_distances = np.zeros(num_samples, dtype=float)
    success_bools = np.zeros(num_samples, dtype=bool)
    # N.B. 0 is the same as False if you specify that dtype is boolean

    # For each sample, perform optimization
    for sample_ind in range(num_samples):
        sample_vals = data_df.iloc[:, sample_ind].values

        # Calculate optimized distances
        optimization_result = minimize_scalar(distance_function, args=(sample_vals, probe_medians),
                                              method="bounded", bounds=optim_bounds)

        # Return offset, optimized distance, and a flag indicating if the
        # solution converged
        optimized_offsets[sample_ind] = optimization_result.x
        optimized_distances[sample_ind] = optimization_result.fun
        success_bools[sample_ind] = optimization_result.success

    # Apply the optimized offsets
    df_with_offsets = data_df.add(optimized_offsets)

    # Report which samples did not converge
    if not all(success_bools):
        logger.info(("The following samples failed to converge during " +
                     "optimization: \n{}").format(data_df.columns[~success_bools].values))

    # Return the df with offsets applied, the offsets, the optimized distances,
    # and boolean array of samples that converged
    return df_with_offsets, optimized_offsets, optimized_distances, success_bools


# tested #
def calculate_distances_only(data_df):
    """For each sample, calculate distance metric.
    N.B. Only uses the non-NaN values from each sample.

    Args:
        data_df (pandas df)
    Returns:
        dists (numpy array of floats): length = num_samples
    """
    # Determine the median value for each probe
    probe_medians = data_df.median(axis=1)

    # Initialize output
    num_samples = data_df.shape[1]
    dists = np.zeros(num_samples, dtype=float)

    # For each sample, calculate distance
    for sample_ind in range(num_samples):
        sample_vals = data_df.iloc[:, sample_ind].values

        # Calculate unoptimized distances
        dists[sample_ind] = distance_function(0, sample_vals, probe_medians)

    return dists


# tested #
def distance_function(offset, values, medians):
    """This function calculates the distance metric.
    N.B. Only uses the non-NaN values.

    Args:
        offset (float)
        values (numpy array of floats)
        medians (numpy array of floats)
    Returns:
        dist (float)
    """
    non_nan_idx = ~np.isnan(values)
    assert np.size(non_nan_idx) != 0, "All values in this sample are NaN!"

    non_nan_values = values[non_nan_idx]
    non_nan_medians = medians[non_nan_idx]
    dist = sum(np.square((offset + non_nan_values) - non_nan_medians))
    return dist


def p100_filter_samples_by_dist(gct, assay_type, p100_assay_types, offsets,
                                dists, success_bools, dist_sd_cutoff, prov_code):
    """If P100, filter out samples whose distance metric is above some threshold.
    Also remove samples for which optimization did not converge.

    N.B. Offsets are only passed to this function so that offsets of filtered
    out samples can be removed.

    Args:
        gct (GCToo object)
        assay_type (string)
        p100_assay_types (list of strings)
        offsets (numpy array of floats, or None): offset that was applied to each sample
        dists (numpy array of floats, or None): distance metric for each sample
        success_bools (numpy array of bools, or None): for each sample, indicates whether optimization was successful
        dist_sd_cutoff (float): maximum SD for a sample's distance metric before being filtered out
        prov_code (list of strings)

    Returns:
        gct (GCToo object): updated
        out_offsets (numpy array of floats, or None): only returned for samples that were not removed
        prov_code (list of strings): updated

    """
    # P100
    if (assay_type in p100_assay_types):
        (out_df, out_offsets) = remove_sample_outliers(gct.data_df, offsets, dists, success_bools, dist_sd_cutoff)
        prov_code_entry = OUTLIER_SAMPLE_FILTER_PROV_CODE_ENTRY
        # prov_code_entry = "OF{:.0f}".format(args.dist_sd_cutoff) # TODO(lev): should we include dist_sd_cutoff in this prov_code_entry?
        (gct, prov_code) = update_metadata_and_prov_code(
            out_df, gct.row_metadata_df, gct.col_metadata_df, prov_code_entry, prov_code)

    # GCP
    else:
        out_offsets = None

    return gct, out_offsets, prov_code


# tested #
def remove_sample_outliers(data_df, offsets, distances, success_bools, dist_sd_cutoff):
    """Calculate distance cutoff for outlier samples and remove outliers.

    Samples are considered outliers if their distance is above the cutoff
    OR if the optimization process didn't converge for them. Return only the
    offsets of samples that were not removed.

    Args:
        data_df (pandas df)
        offsets (numpy array of floats): length = num_samples
        distances (numpy array of floats): length = num_samples
        success_bools (numpy array of bools): length = num_samples
        dist_sd_cutoff (float): maximum SD for a sample's distance metric before being filtered out

    Returns:
        out_df (pandas df): potentially smaller than original df
        out_offsets (numpy array of floats): length = num_remaining_samples
    """
    assert len(distances) == data_df.shape[1], (
        "len(distances): {} does not equal data_df.shape[1]: {}").format(len(distances), data_df.shape[1])

    # Calculate the distance cutoff
    cutoff = np.mean(distances) + np.multiply(np.std(distances), dist_sd_cutoff)

    # If optimization occurred, consider success_bools; otherwise, ignore them
    if success_bools is not None:
        # Keep only samples whose distance metric is less than the cutoff AND
        # converged during optimization
        samples_to_keep = np.logical_and(distances < cutoff, success_bools)
    else:
        samples_to_keep = distances < cutoff

    out_df = data_df.iloc[:, samples_to_keep]
    assert not out_df.empty, "All samples were filtered out. Try increasing the SD cutoff."

    # If optimization occurred, return only subsets of remaining samples
    if offsets is not None:
        out_offsets = offsets[samples_to_keep]
    else:
        # Otherwise, return None for the offsets
        out_offsets = None

    return out_df, out_offsets


def median_normalize(gct, subset_normalize_bool, row_subset_field, col_subset_field, prov_code):
    """Subset normalize if subset_normalize_bool=True OR if the metadata
    shows that subsets exist for either rows or columns. Otherwise, use the
    whole row for median normalization.

    Args:
        gct (GCToo object)
        subset_normalize_bool (bool): true indicates that subset normalization should be performed
        row_subset_field (string): row metadata field indicating the row subset group
        col_subset_field (string): col metadata field indicating the col subset group
        prov_code (list of strings)

    Returns:
        gct (GCToo object): updated
        updated_prov_code (list of strings)

    """
    # Check if subsets_exist
    subsets_exist = check_for_subsets(gct.row_metadata_df, gct.col_metadata_df,
                                      row_subset_field, col_subset_field)

    # TODO(lev): should this be AND or OR?
    if (subset_normalize_bool) and (subsets_exist):
        logger.debug("Performing subset normalization.")
        out_gct = subset_normalize(gct, row_subset_field, col_subset_field)
        prov_code_entry = SUBSET_NORMALIZE_PROV_CODE_ENTRY
        updated_prov_code = prov_code + [prov_code_entry]
    else:
        # Subtract median of whole row from each entry in the row
        gct.data_df = row_median_normalize(gct.data_df)
        out_gct = gct
        prov_code_entry = ROW_NORMALIZE_PROV_CODE_ENTRY
        updated_prov_code = prov_code + [prov_code_entry]

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
def subset_normalize(gct, row_subset_field, col_subset_field):
    """Row-median normalize subsets of data.

    The function make_norm_ndarray extracts relevant metadata to figure out
    which subsets of the data should be row median normalized together, and
    iterate_over_norm_ndarray_and_normalize actually applies the normalization.

    Args:
        gct (GCToo object)
        row_subset_field (string): row metadata field indicating the row subset group
        col_subset_field (string): col metadata field indicating the col subset group

    Returns:
        gct (GCToo object): with normalized data
    """
    # Create normalization ndarray
    norm_ndarray = make_norm_ndarray(gct.row_metadata_df, gct.col_metadata_df, row_subset_field, col_subset_field)

    # Iterate over norm ndarray and actually perform normalization
    gct.data_df = iterate_over_norm_ndarray_and_normalize(gct.data_df, norm_ndarray)

    return gct

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

    # Convert sample group vectors from strs to lists of strings
    sample_grps_lists = [sample_str.split(",") for sample_str in sample_grps_strs]

    # Verify that all lists have same length
    length_of_first_list = len(sample_grps_lists[0])
    assert all([len(sample_list) == length_of_first_list for sample_list in sample_grps_lists])

    # Convert from lists of strings to ndarray of ints
    sample_grp_ndarray = np.array(sample_grps_lists, dtype="int")

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
    assert ~np.any(norm_ndarray == 0)

    return norm_ndarray

# tested #
def iterate_over_norm_ndarray_and_normalize(data_df, norm_ndarray):
    """Iterate over norm_ndarray and row-median normalize the subsets.

    Each row of norm_ndarray indicates the subsets for normalization to
    select from data_df.

    Args:
        data_df (pandas df): size = (num_probes, num_samples)
        norm_ndarray (numpy ndarray): size = (num_probes, num_samples)
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

            # Subtract median and insert into output ndarray
            data_to_insert = data_to_normalize - np.nanmedian(data_to_normalize)
            normalized_data.iloc[row_idx].loc[norm_ndarray_row == sample_subset] = data_to_insert

    return normalized_data

# tested #
def row_median_normalize(data_df):
    """Subtract median of the row from each entry of the row.

    Args:
        data_df (pandas df)
    Returns:
        out_df (pandas df): with normalized values
    """
    out_df = data_df.subtract(data_df.median(axis=1), axis=0)
    return out_df


# tested #
def insert_offsets_and_prov_code(gct, offsets, offsets_field, prov_code, prov_code_field, prov_code_delimiter):
    """Insert offsets into output gct and update provenance code in metadata.

    Args:
        gct (GCToo object)
        offsets (numpy array of floats, or None): if optimization was not performed, this will be None
        offsets_field (string): name of col metadata field into which offsets will be inserted
        prov_code (list of strings)
        prov_code_field (string): name of col metadata field containing the provenance code
        prov_code_delimiter (string): what string to use as delimiter in prov_code

    Returns:
        gct (GCToo object): updated metadata

    """
    # Make sure that metadata and data dfs agree in shape
    assert gct.data_df.shape[0] == gct.row_metadata_df.shape[0]
    assert gct.data_df.shape[1] == gct.col_metadata_df.shape[0]

    # If offsets is not none, insert into col_metadata_df
    if offsets is not None:
        assert len(offsets) == gct.col_metadata_df.shape[0]
        gct.col_metadata_df[offsets_field] = offsets.astype(str)

    # TODO(lev): check that only acceptable values were inserted?

    # Convert provenance code to delimiter separated string
    prov_code_str = prov_code_delimiter.join(prov_code)

    # Update the provenance code in col_metadata_df
    gct.col_metadata_df.loc[:, prov_code_field] = prov_code_str

    return gct


def write_output_gct(gct, out_path, out_name, data_null, filler_null):
    """Write output gct file.

    Args:
        gct (GCToo object)
        out_path (string): path to save directory
        out_name (string): name of output gct
        data_null (string): string with which to represent NaN in data
        filler_null (string): string with which to fill the empty top-left quadrant in the output gct

    Returns:
        None

    """
    out_fname = os.path.join(out_path, out_name)
    write_gctoo.write(gct, out_fname, data_null=data_null,
                      filler_null=filler_null, data_float_format=None)


def update_metadata_and_prov_code(data_df, row_meta_df, col_meta_df, prov_code_entry, prov_code):
    """Update metadata with the already sliced data_df, and update the prov_code.

    Args:
        data_df (pandas df)
        row_meta_df (pandas df)
        col_meta_df (pandas df)
        prov_code_entry (string)
        prov_code (list_of_strings)

    Returns:
        out_gct (GCToo object): updated
        updated_prov_code (list of strings)

    """
    if prov_code_entry is not None:
        updated_prov_code = prov_code + [prov_code_entry]
    else:
        updated_prov_code = prov_code

    out_gct = slice_metadata_using_already_sliced_data_df(data_df, row_meta_df, col_meta_df)
    return out_gct, updated_prov_code


# tested #
def slice_metadata_using_already_sliced_data_df(data_df, row_meta_df, col_meta_df):
    """Slice row_meta_df and col_meta_df to only contain the row_ids and col_ids
    in data_df.

    Args:
        data_df (pandas df)
        row_meta_df (pandas df)
        col_meta_df (pandas df)

    Returns:
        out_gct (GCToo object): with all dfs correctly sliced

    """
    # Get rows and samples to keep from data_df
    rows = data_df.index.values
    cols = data_df.columns.values

    # Number of rows and samples should be > 1
    # If only one row or col, the df will be transposed and that will be ugly
    if not len(rows) > 1:
        err_msg = "Fewer than 2 rows remain after data processing. I don't like that!"
        logger.error(err_msg)
        raise Exception(err_msg)
    if not len(cols) > 1:
        err_msg = "Fewer than 2 columns remain after data processing. I don't like that!"
        logger.error(err_msg)
        raise Exception(err_msg)

    # Extract remaining rows and samples from row and col metadata dfs
    out_row_meta_df = row_meta_df.loc[rows, :]
    out_col_meta_df = col_meta_df.loc[cols, :]

    # Return gct
    out_gct = GCToo.GCToo(data_df=data_df,
                          row_metadata_df=out_row_meta_df,
                          col_metadata_df=out_col_meta_df)

    return out_gct


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args: {}".format(args))

    # Check that config file exists
    assert os.path.exists(os.path.expanduser(args.psp_config_path))

    main(args)
