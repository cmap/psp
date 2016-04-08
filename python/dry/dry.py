import logging
import utils.setup_logger as setup_logger
import ConfigParser
import argparse
import sys
import os

import numpy as np
from scipy.optimize import minimize_scalar
import pandas as pd
import in_out.GCToo as GCToo
import in_out.parse_gctoo as parse_gctoo
import in_out.write_gctoo as write_gctoo

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

""" Performs filtering and normalization of p100 and GCP data.

Converts level 2 to level 3 data. Required input is a filepath to a gct file.
Output is a gct file.

N.B. This script requires a configuration file in your home directory
(e.g. ~/.PSP_config).

"""

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required argument
    parser.add_argument("gct_path", help="filepath to gct file", type=str)

    # Optional arguments
    parser.add_argument("-out_path", default=None, type=str,
                        help="directory where to save the output files")
    parser.add_argument("-out_name", default=None, type=str,
                        help="basename of the output files")
    parser.add_argument("-PSP_config_path", type=str,
                        default="~/.PSP_config",
                        help="path to PSP config file")
    parser.add_argument("-process_mode", choices=["full", "quick"], # doesn't do anything yet
                        default="full",
                        help="whether to run a quick version of the process")
    parser.add_argument("-log2", action="store_true",
                        help="whether to perform log2 normalization")

    # TODO(lev): defaults should be different depending on the assay. How to do this...?

    parser.add_argument("-sample_nan_thresh", type=float, default=0.3,
                        help=("if < sample_nan_thresh of a sample's data " +
                              "is non-nan, that sample will be filtered out"))
    parser.add_argument("-probe_nan_thresh", type=float, default=0.3,
                        help=("if < probe_nan_thresh of a probe's data " +
                              "is non-nan, that probe will be filtered out"))
    parser.add_argument("-probe_sd_cutoff", type=float, default=3,
                        help=("maximum SD for a probe " +
                              "before being filtered out"))
    parser.add_argument("-dist_sd_cutoff", type=float, default=3,
                        help=("maximum SD for a sample's distance metric " +
                              "before being filtered out"))
    parser.add_argument("-prov_code_delimiter", type=str, default="+",
                        help=("string delimiter in provenance code"))
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                        help="increase the number of messages reported")

    # P100-specific arguments
    parser.add_argument("-optim", action="store_true",
                        help="whether to perform load balancing optimization")
    parser.add_argument("-optim_bounds", type=tuple, default=(-7,7),
                        help="bounds over which to perform optimization")

    # GCP-specific arguments
    parser.add_argument("-normalization_peptide_id", type=str, default="BI10052",
                        help=("which peptide to normalize to"))
    parser.add_argument("-subset_normalize_bool", action="store_true",
                        help="whether to perform probe group normalization")

    # TODO(lev): add argument "force_assay"? (JJ's suggestion)

    return parser


def main(args):
    """THE MAIN METHOD. Filters and normalizes data and save result as a gct file.

    Args:
        args (argparse.Namespace): fields as defined in build_parser()
    Returns:
        out_gct (GCToo): output gct object
    """
    # 1) Read gct and extract provenance code
    (gct, prov_code) = read_gct_and_check_provenance_code(args.PSP_config_path, args.gct_path)

    # 2) Log transform and, if GCP, histone normalize data
    # (gct.data_df, prov_code) = log_and_gcp_histone_normalize(
    #     gct.data_df, args.normalization_peptide_id, prov_code)

    ### LOG TRANSFORM
    (log_transformed_df, prov_code) = do_log_transform_if_needed(gct.data_df, prov_code)

    ### GCP HISTONE NORMALIZE (if GCP)
    (out_df, prov_code) = do_gcp_histone_normalize_if_needed(log_transformed_df, args.normalization_peptide_id, prov_code)

    ##########

    # 3) Filter out samples
    # (gct.data_df, prov_code) = filter_samples(
    #     gct.data_df, args.sample_nan_thresh, args.optim, args.optim_bounds,
    #     args.dist_sd_cutoff, prov_code)

    ### FILTER SAMPLES BY NAN
    gct.data_df = filter_samples_by_nan(gct.data_df, args.sample_nan_thresh)
    thresh_digit = ("{:.1f}".format(args.sample_nan_thresh)).split(".")[1]
    prov_code_entry = "SF{}".format(thresh_digit)
    prov_code.append(prov_code_entry)

    # 4) Filter out probes
    # (gct.data_df, prov_code) = filter_probes(
    #     gct.data_df, gct.row_metadata_df, prov_code, args.probe_nan_thresh,
    #     args.probe_sd_cutoff)

    ### FILTER MANUALLY REJECTED PROBES
    gct.data_df = manual_probe_rejection(gct.data_df, gct.row_metadata_df)
    # prov_code_entry = "MPR"
    # prov_code.append(prov_code_entry)

    # TODO(lev): I added this provenance code. Do you want to keep it?

    ### FILTER PROBES BY NAN AND SD
    gct.data_df = filter_probes_by_nan_and_sd(gct.data_df, args.probe_nan_thresh, args.probe_sd_cutoff)
    thresh_digit = ("{:.1f}".format(args.probe_nan_thresh)).split(".")[1]
    prov_code_entry = "PF{}".format(thresh_digit)
    prov_code.append(prov_code_entry)

    ##########

    # TODO(lev): should only happen for P100
    ### CALCULATE DISTANCES
    (gct.data_df, offsets, dists, success_bools, prov_code) = calculate_distances_and_optimize_if_needed(gct.data_df, args.optim, args.optim_bounds, prov_code)

    ### REMOVE SAMPLE OUTLIERS
    gct.data_df = remove_sample_outliers(gct.data_df, dists, success_bools, args.dist_sd_cutoff)
    prov_code_entry = "OSF"
    prov_code.append(prov_code_entry)

    ##########

    # 5) Row median and probe-group specific normalize
    # (gct.data_df, prov_code) = row_median_and_subset_normalize(gct.data_df, prov_code, args.subset_normalize)

    ### ROW MEDIAN NORMALIZE
    gct.data_df = row_median_normalize(gct.data_df)
    prov_code_entry = "RMN"
    prov_code.append(prov_code_entry)

    # TODO(lev): make sure this works.

    ### SUBSET NORMALIZE
    if args.subset_normalize_bool:
        (gct.data_df, prov_code) = subset_normalize(gct.data_df, gct.row_metadata_df, gct.col_metadata_df, prov_code)

    ##########

    # 6) Produce QC figures
    # Does Plato handle non-384 well plates?
    # One pixel is value, one pixel is error-bar?

    # 7) Create processed gct object
    out_gct = create_output_gct(gct.data_df, gct.row_metadata_df, gct.col_metadata_df,
                                offsets, prov_code, args.prov_code_delimiter)

    # 8) Save processed gct object to file
    out_fname = os.path.join(args.out_path, args.out_name)
    write_gctoo.write(out_gct, out_fname, data_null="NaN", filler_null="NA", data_float_format=None)
    return out_gct


# tested #
def read_gct_and_check_provenance_code(PSP_config_path, gct_path):
    """ Create gct object and check that provenance code is non-empty and the same for all samples."""

    # Read from config file
    configParser = ConfigParser.RawConfigParser()
    configParser.read(os.path.expanduser(PSP_config_path))

    # Extract what values to consider NaN
    # N.B. eval used to convert string from config file to list
    PSP_nan_values = eval(configParser.get("io", "nan_values"))

    # Parse the gct file and return gct object
    gct = parse_gctoo.parse(gct_path, nan_values=PSP_nan_values)

    # Extract the plate's provenance code and assay type
    prov_code = extract_prov_code(gct.col_metadata_df)
    assay_type = prov_code[0]

    # TODO(lev): check that prov code only has allowed values

    # Make sure assay_type is one of the allowed values
    P100_assay_types = eval(configParser.get("metadata", "P100_assays"))
    GCP_assay_types = eval(configParser.get("metadata", "GCP_assays"))
    allowed_assay_types = P100_assay_types + GCP_assay_types
    check_assay_type(assay_type, allowed_assay_types)

    return gct, prov_code


def log_and_gcp_histone_normalize(data_df, normalization_peptide_id, prov_code):
    """Perform log transformation and, if GCP assay, histone normalization.

    Args:
        data_df: dataframe of floats
        normalization_peptide_id: string of peptide to normalize to
        prov_code: list of strings
    Returns:
        out_df: dataframe of floats with one row removed if histone normalized
        prov_code: updated provenance code
    """
    ### LOG TRANSFORM
    (log_transformed_df, prov_code) = do_log_transform_if_needed(data_df, prov_code)

    ### GCP HISTONE NORMALIZE (if GCP)
    (out_df, prov_code) = do_gcp_histone_normalize_if_needed(
        log_transformed_df, normalization_peptide_id, prov_code)

    return out_df, prov_code


def filter_samples(data_df, sample_nan_thresh, optim, optim_bounds, dist_sd_cutoff, prov_code):

    ### FILTER SAMPLES BY NAN
    data_df = filter_samples_by_nan(data_df, sample_nan_thresh)
    thresh_digit = ("{:.1f}".format(sample_nan_thresh)).split(".")[1]
    prov_code_entry = "SF{}".format(thresh_digit)
    prov_code.append(prov_code_entry)

    # TODO(lev): should only happen for P100?
    ### CALCULATE DISTANCES
    (data_df, dists, success_bools, prov_code) = (
        calculate_distances_and_optimize_if_needed(
            data_df, optim, optim_bounds, prov_code))

    # TODO(lev): should only happen for P100?
    ### REMOVE SAMPLE OUTLIERS
    data_df = remove_sample_outliers(data_df, dists, success_bools, dist_sd_cutoff)
    prov_code_entry = "OSF"
    prov_code.append(prov_code_entry)

    return data_df, prov_code

# probably need to separate this into diff fcns
def filter_probes(data_df, row_metadata_df, prov_code, probe_nan_thresh, probe_sd_cutoff):

    ### FILTER MANUALLY REJECTED PROBES
    data_df = manual_probe_rejection(data_df, row_metadata_df)
    prov_code_entry = "MPR"
    prov_code.append(prov_code_entry)

    # TODO(lev): I added this provenance code. Do you want to keep it?

    ### FILTER PROBES BY NAN AND SD
    data_df = filter_probes_by_nan_and_sd(data_df, probe_nan_thresh, probe_sd_cutoff)
    thresh_digit = ("{:.1f}".format(probe_nan_thresh)).split(".")[1]
    prov_code_entry = "PF{}".format(thresh_digit)
    prov_code.append(prov_code_entry)

    return data_df, prov_code


def row_median_and_subset_normalize(data_df, prov_code, subset_normalize_bool):

    ### ROW MEDIAN NORMALIZE
    data_df = row_median_normalize(data_df)
    prov_code_entry = "RMN"
    prov_code.append(prov_code_entry)

    ### SUBSET NORMALIZE
    if subset_normalize_bool:
        (data_df, prov_code) = subset_normalize(data_df,row_metadata_df, col_metadata_df, prov_code)

    return data_df, prov_code


def create_output_gct(data_df, row_df, col_df, offsets, prov_code, prov_code_delimiter):
    """Create gct object from component dataframes.

    Uses probe and sample ids from data_df to return correct row and col metadata.
    Also, reinserts the updated provenance code and appends correction_factor as
    a column metadata header.

    Args:
        data_df (pandas df)
        row_df  (pandas df)
        col_df (pandas df)
        offsets (numpy array or None): if optim, then type will be numpy array; otherwise None
        prov_code (list of strings)
        prov_code_delimiter (string): what string to use as delimiter in prov_code
    Returns:
        out_gct (GCToo): output gct
    """
    PROV_CODE_FIELD = "provenance_code"

    # TODO(lev): clean this up. How do I make sure that len(offsets) == data_df.shape[1] ??

    # Insert offsets as string field in col_metadata_df unless it's None
    if offsets is not None:
        col_df["optimization_offset"] = offsets.astype(str)

    # Get remaining rows and samples from data_df
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
    out_row_df = row_df.loc[rows, :]
    out_col_df = col_df.loc[cols, :]

    # Verify shapes of output dfs
    assert data_df.shape[0] == out_row_df.shape[0]
    assert data_df.shape[1] == out_col_df.shape[0]

    # TODO(lev): check that only acceptable values were inserted

    # Convert provenance code to delimiter separated string
    prov_code_str = prov_code_delimiter.join(prov_code)

    # Update the provenance code in col_metadata_df
    out_col_df.loc[:, PROV_CODE_FIELD] = prov_code_str

        # Insert component dfs into gct object
    out_gct = GCToo.GCToo(row_metadata_df=out_row_df,
                          col_metadata_df=out_col_df,
                          data_df=data_df)
    return out_gct


# tested #
def extract_prov_code(col_metadata_df):
    """Extract the provenance code from the column metadata.
    Also verify that it is the same for all samples.

    Args:
        col_metadata_df (pandas df): contains provenance code metadata
    Returns:
        prov_code (list of strings)
    """
    # Create pandas series of all provenance codes
    prov_code_series = col_metadata_df.loc[:, "provenance_code"]

    # Split each provenance code string along the separator
    prov_code_list_series = prov_code_series.apply(lambda x: x.split("+"))

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


def do_log_transform_if_needed(data_df, prov_code):
    """Call log_transform if it hasn't already been done."""

    # Check if log2 transformation has already occurred
    if "L2X" in prov_code:
        logger.info("L2X has already occurred.")
        return data_df, prov_code
    else:
        # Perform log2 transformation and update provenance code
        out_df = log_transform(data_df, log_base=2)
        prov_code_entry = "L2X"
        prov_code.append(prov_code_entry)
        return out_df, prov_code

# tested #
def log_transform(data_df, log_base):
    """Take the log of each value in a pandas dataframe.

    Args:
        data_df: pandas dataframe of floats
        log_base: the base of the logarithm
    Returns:
        out_df: pandas dataframe with log transformed values
    """
    # Replace 0 with np.nan
    data_df.replace(0, np.nan, inplace=True)

    # Numpy operations work fine on dataframes
    out_df = np.log(data_df) / np.log(log_base)
    return out_df


def do_gcp_histone_normalize_if_needed(data_df, normalization_peptide_id, prov_code):
    """If GCP assay, call gcp_histone_normalize. """

    assay_type = prov_code[0]


    # TODO(lev): move this to config file too
    if assay_type in ["GR1"]:
        data_df = gcp_histone_normalize(data_df, normalization_peptide_id)
        prov_code_entry = "H3N"
        prov_code.append(prov_code_entry)
    return data_df, prov_code

# tested #
def gcp_histone_normalize(data_df, normalization_peptide_id):
    """Subtract values of normalization_peptide_id from all values.
    Also, remove the row of data corresponding to the normalization histone.

    Args:
        data_df: dataframe of floats
        normalization_peptide_id: string
    Returns:
        out_df: dataframe of floats with one row removed
    """
    # Verify that normalization peptide is in the data
    assert normalization_peptide_id in data_df.index, (
        ("The normalization peptide is not in this dataset. " +
         "normalization_peptide_id: {}".format(normalization_peptide_id)))

    # Calculate normalization values
    norm_values = data_df.loc[normalization_peptide_id, :]

    # Drop the normalization peptide row
    data_df.drop(normalization_peptide_id, inplace=True)

    # Subtract the normalization values from all rows
    out_df = data_df - norm_values
    return out_df

# tested #
def filter_samples_by_nan(data_df, sample_nan_thresh):
    """Remove samples (i.e. columns) with less than sample_nan_thresh non-NaN values.

    Args:
        data_df: pandas dataframe of floats
        sample_nan_thresh: float from 0 to 1
    Returns:
        out_df: pandas dataframe (potentially smaller than original df)
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
def manual_probe_rejection(data_df, row_metadata_df):
    """Remove probes that were manually selected for rejection.

    Args:
        data_df: pandas dataframe of floats
        row_metadata_df: pandas dataframe of strings
    Returns:
        out_df: pandas dataframe of floats (potentially smaller than input)
    """
    # Extract "pr_probe_suitability_manual" metadata field
    keep_probe_str = row_metadata_df.loc[:, "pr_probe_suitability_manual"]

    # Convert strings to booleans
    keep_probe_bool = (keep_probe_str == "TRUE")

    # Check that the list of good probes is not empty
    assert keep_probe_bool.any(), ("No probes were marked TRUE (i.e. suitable).\n" +
                                   "row_metadata_df.loc[:, " +
                                   "'pr_probe_suitability_manual']: \n{}").format(keep_probe_str)

    # Return the good probes
    out_df = data_df[keep_probe_bool.values]
    assert not out_df.empty, "All probes were filtered out!"
    return out_df

# tested #
def filter_probes_by_nan_and_sd(data_df, probe_nan_thresh, probe_sd_cutoff):
    """Remove probes (i.e. rows) with less than probe_nan_thresh non-NaN values.
    Also remove probes with standard deviation higher than probe_sd_cutoff.

    Args:
        data_df: pandas dataframe of floats
        probe_nan_thresh: float from 0 to 1
        probe_sd_cutoff: float
    Returns:
        out_df: pandas dataframe (potentially smaller than original df)
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
    assert not out_df.empty, ("All probes were filtered out. " +
                              "Try reducing the NaN threshold and/or SD cutoff.")
    return out_df

# tested #
def calculate_distances_and_optimize_if_needed(data_df, optim_bool, optim_bounds, prov_code):
    """If P100 assay and optim=True, call calculate_distances_and_optimize.
    Otherwise, call calculate_distances.
    """
    assay_type = prov_code[0]

    # TODO(lev): P100_assays should go into config file

    if assay_type in ["PR1", "DIA1", "PRM"] and optim_bool:
        (data_df, offsets, dists, success_bools) = calculate_distances_and_optimize(data_df, optim_bounds)
        prov_code_entry = "LLB"
        prov_code.append(prov_code_entry)
    else:
        # Simply calculate distance metric for each sample
        dists = calculate_distances_only(data_df)

        # Create artificial success_bools to be used by remove_sample_outliers
        success_bools = np.ones(len(dists), dtype=bool)

        # TODO(lev): This should be better. Simply return None for offsets
        offsets = None

    return data_df, offsets, dists, success_bools, prov_code

# tested #
def calculate_distances_and_optimize(data_df, optim_bounds):
    """For each sample, perform optimized load balancing.

    This means finding an offset for each sample that minimizes
    the distance metric. The distance metric is the sum of the
    distances from each probe measurement to its median.

    N.B. Only uses the non-NaN values from each sample.

    Args:
        data_df: pandas dataframe of floats
        optim_bounds: tuple of floats

    Returns:
        out_df: pandas dataframe with offsets applied
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
        data_df: pandas dataframe of floats
    Returns:
        dists: numpy array with length = num_samples
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
        offset: float
        values: numpy array of floats
        medians: numpy array of floats
    Returns:
        dist: float
    """
    non_nan_idx = ~np.isnan(values)
    assert np.size(non_nan_idx) != 0, "All values in this sample are NaN!"

    non_nan_values = values[non_nan_idx]
    non_nan_medians = medians[non_nan_idx]
    dist = sum(np.square((offset + non_nan_values) - non_nan_medians))
    return dist

# tested #
def remove_sample_outliers(data_df, distances, success_bools, dist_sd_cutoff):
    """Calculate distance cutoff for outlier samples and remove outliers.

    Samples are considered outliers if their distance is above the cutoff
    OR if the optimization process didn't converge for them.

    Args:
        data_df: pandas dataframe of floats
        distances: numpy array of floats with length = num_samples
        success_bools: numpy array of bools with length = num_samples
        dist_sd_cutoff: float
    Returns:
        out_df: pandas dataframe (potentially smaller than original df)
    """
    assert len(distances) == data_df.shape[1], (
        "len(distances): {} does not equal " +
        "data_df.shape[1]: {}").format(len(distances), data_df.shape[1])
    assert len(success_bools) == data_df.shape[1], (
        "len(success_bools): {} does not equal " +
        "data_df.shape[1]: {}").format(len(success_bools), data_df.shape[1])

    # Calculate the distance cutoff
    cutoff = np.mean(distances) + np.multiply(np.std(distances),
                                              dist_sd_cutoff)

    # Remove samples whose distance metric is greater than the cutoff OR
    # those that didn't converge during optimization
    out_df = data_df.iloc[:, np.logical_and(distances < cutoff, success_bools)]
    assert not out_df.empty, "All samples were filtered out. Try reducing the SD cutoff."
    return out_df

# tested #
def row_median_normalize(data_df):
    """Subtract the row median from values.

    Args:
        data_df: dataframe of floats
    Returns:
        out_df: dataframe of floats
    """
    out_df = data_df.subtract(data_df.median(axis=1), axis=0)
    return out_df

# tested #
def subset_normalize(data_df, row_metadata_df, col_metadata_df, prov_code):
    """Row-median normalize subsets of data.

    The function make_norm_ndarray extracts relevant metadata to figure out
    which subsets of the data should be row median normalized together, and
    iterate_over_norm_ndarray_and_normalize actually applies the normalization.

    Args:
        data_df (pandas df)
        row_metadata_df (pandas df)
        col_metadata_df (pandas df)
        prov_code (list of strings)
    Returns:
        data_df (pandas df): with normalized data
        prov_code (list of strings): updated provenance code
    """
    # Create normalization ndarray
    norm_ndarray = make_norm_ndarray(row_metadata_df, col_metadata_df)

    # Iterate over norm ndarray and actually perform normalization
    data_df = iterate_over_norm_ndarray_and_normalize(data_df, norm_ndarray)

    # Update provenance code
    prov_code_entry = "GMN"
    prov_code.append(prov_code_entry)

    return data_df, prov_code

# tested #
def make_norm_ndarray(row_metadata_df, col_metadata_df):
    """Creates a normalization ndarray.

    The normalization ndarray is used to figure out how to normalize each
    section of the data. Extracts probe normalization groups from row
    metadata field "pr_probe_normalization_group" and sample normalization
    groups from column metadata field "det_normalization_group_vector."

    Args:
        row_metadata_df (pandas df)
        col_metadata_df (pandas df)
    Returns:
        norm_ndarray (numpy array): size = (num_probes, num_samples)
    """
    PROBE_GROUP_FIELD = "pr_probe_normalization_group"
    SAMPLE_GROUP_FIELD = "det_normalization_group_vector"

    # Verfiy that metadata fields of interest exist
    assert PROBE_GROUP_FIELD in row_metadata_df.columns
    assert SAMPLE_GROUP_FIELD in col_metadata_df.columns

    # Get sample group vectors
    sample_grps_strs = col_metadata_df[SAMPLE_GROUP_FIELD].values

    # Convert sample group vectors from strs to lists of strings
    sample_grps_lists = [sample_str.split(",") for sample_str in sample_grps_strs]

    # Verify that all lists have same length
    length_of_first_list = len(sample_grps_lists[0])
    assert all([len(sample_list) == length_of_first_list for sample_list in sample_grps_lists])

    # Convert from lists of strings to ndarray of ints
    sample_grp_ndarray = np.array(sample_grps_lists, dtype='int')

    # Get probe groups and unique probe groups; convert to ints
    probe_grps = row_metadata_df[PROBE_GROUP_FIELD].values.astype('int')
    unique_probe_grps = np.unique(probe_grps)

    # Initialize norm_ndarray
    norm_ndarray = np.zeros((row_metadata_df.shape[0], col_metadata_df.shape[0]), dtype='int')

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


if __name__ == '__main__':
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args: {}".format(args))

    # Check that config filepath exists
    assert os.path.exists(os.path.expanduser(args.PSP_config_path))

    main(args)
