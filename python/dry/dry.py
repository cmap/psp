import logging
import utils.setup_logger as setup_logger
import ConfigParser
import argparse
import sys

import numpy as np
from scipy.optimize import minimize_scalar
import pandas as pd
import in_out.parse_gctoo as parse_gctoo

logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Read config file
PSP_config_path = "/Users/lev/.PSP_config"
configParser = ConfigParser.RawConfigParser()
configParser.read(PSP_config_path)


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required argument
    parser.add_argument("gct_path", help="filepath to gct file", type=str)

    # Optional arguments
    parser.add_argument("-out_path", default=None,
                        help="directory where to save the output files")
    parser.add_argument("-out_name", default=None,
                        help="basename of the output files")
    parser.add_argument("-process_mode", choices=["full", "quick"],
                        default="full",
                        help="whether to run a quick version of the process")
    parser.add_argument("-log2", action="store_true",
                        help="whether to perform log2 normalization")
    parser.add_argument("-sample_nan_thresh", type=float, default=0.3,
                        help=("if < sample_nan_thresh of a sample's data " +
                              "is non-nan, that sample will be filtered out"))
    parser.add_argument("-sample_nan_thresh", type=float, default=0.3,
                        help=("if < probe_nan_thresh of a probe's data " +
                              "is non-nan, that probe will be filtered out"))
    parser.add_argument("-probe_sd_cutoff", type=float, default=3,
                        help=("maximum SD for a probe " +
                              "before being filtered out"))
    parser.add_argument("-dist_sd_cutoff", type=float, default=3,
                        help=("maximum SD for a sample's distance metric " +
                              "before being filtered out"))
    parser.add_argument("-run", action="store_true", default=False,
                        help="whether to actually run the command")
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
    parser.add_argument("-probe_group_normalize", action="store_true",
                        help="whether to perform probe group normalization")

    return parser


def main(args):
    """The main method. Writes processed gct to a file and produces QC figures.

    Args:
        args: object with fields as defined in build_parser()
    Returns:
        out_gct: gct object post-processing
    """
    # Extract what values to consider NaN from config file
    # N.B. eval used to convert string from config file to list
    PSP_nan_values = eval(configParser.get("io", "nan_values"))

    # Parse the gct file and return gct object
    gct = parse_gctoo.parse(args.gct_path, nan_values=PSP_nan_values)

    # Extract the plate's provenance code and assay type
    prov_code = extract_prov_code(gct.col_metadata_df)
    assay_type = prov_code[0]

    # TO-DO(lev): check prov code only has allowed values

    # Make sure assay_type is one of the allowed values
    allowed_assay_types = eval(configParser.get("metadata", "allowed_assay_types"))
    assay_ok = check_assay_type(assay_type, allowed_assay_types)

    ### LOG TRANSFORM
    (gct.data_df, prov_code) = do_log_transform_if_needed(gct.data_df, prov_code)

    ### GCP HISTONE NORMALIZE

    ### FILTER SAMPLES BY NAN
    gct.data_df = filter_samples_by_nan(gct.data_df, args.sample_nan_thresh)
    prov_code_entry = "SF{:.1f}".format(args.sample_nan_thresh).split(".")[1]
    np.append(prov_code, prov_code_entry)

    ### FILTER MANUALLY REJECTED PROBES
    gct.data_df= manual_probe_rejection(gct.data_df, gct.row_metadata_df)
    prov_code_entry = "MPR"
    np.append(prov_code, prov_code_entry)

    ### FILTER PROBES BY NAN AND SD
    gct.data_df = filter_probes_by_nan_and_sd(gct.data_df, args.probe_nan_thresh, args.probe_sd_cutoff)
    prov_code_entry = "PF{:.1f}".format(probe_nan_thresh).split(".")[1]
    np.append(prov_code, prov_code_entry)

    ### CALCULATE DISTANCES
    (data_df, dists, success_bools, prov_code) = (
        calculate_distances_and_optimize_if_needed(
            gct.data_df, assay_type, args.optim, args.optim_bounds, prov_code))

    ### REMOVE SAMPLE OUTLIERS
    gct.data_df = remove_sample_outliers(gct.data_df, dists, success_bools, args.sd_sample_outlier_cutoff)
    prov_code_entry = "OSF"
    np.append(prov_code, prov_code_entry)

    ### ROW MEDIAN NORMALIZE
    gct.data_df = row_median_normalize(gct.data_df)
    prov_code_entry = "RMN"
    np.append(prov_code, prov_code_entry)

    ### PROBE-GROUP SPECIFIC NORMALIZe

    # Save output QC figures

    # Reinsert prov code (check that it only inserted allowed values)
    # Save gct

    # return gct

def extract_prov_code(col_metadata_df):
    """Extract the provenance code from the column metadata.
    Also verify that it is the same for all samples.

    Args:
        col_metadata_df: dataframe containing provenance code metadata
    Returns:
        prov_code: list of strings
    """
    # Create pandas series of all provenance codes
    prov_code_series = col_metadata_df.loc["provenance_code"]

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
        assert prov_code != [""], "Provenance code should not be empty."
        return prov_code
    else:
        err_msg = ("All columns should have the same provenance code, " +
                   "but actually np.unique(prov_code_list_series.values) = {}")
        logger.error(err_msg.format(np.unique(prov_code_list_series.values)))
        raise(Exception(err_msg.format(np.unique(prov_code_list_series.values))))

def check_assay_type(assay_type, allowed_assay_types):
    """Verify that assay is one of the allowed types.

    Args:
        assay_type: string
        allowed_assay_types: list of strings
    Returns:
        assay_ok: bool
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
        out_df = log_transform(data_df)
        prov_code_entry = "L2X"
        np.append(prov_code, prov_code_entry)
        return out_df, prov_code


def log_transform(data_df, log_base=2):
    """Take the log of each value in a pandas dataframe.

    Args:
        data_df: pandas dataframe of floats
        log_base: the base of the logarithm (default = 2).
    Returns:
        out_df: pandas dataframe with log transformed values
    """
    # Replace 0 with np.nan
    data_df.replace(0, np.nan, inplace=True)

    # Numpy operations work fine on dataframes
    out_df = np.log(data_df) / np.log(log_base)
    return out_df


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


def calculate_distances_and_optimize_if_needed(data_df, assay_type, optim_bool, offset_bounds, prov_code):
    """If P100 assay and optim=True, call calculate_distances_and_optimize.
    Otherwise, call calculate_distances.
    """
    if assay_type in ["PR1", "DIA1"] and args.optim:
        (data_df, dists, success_bools) = calculate_distances_and_optimize(data_df, offset_bounds)
        prov_code_entry = "LLB"
        np.append(prov_code, prov_code_entry)
    else:
        # Simply calculate distance metric for each sample
        dists = calculate_distances_only(data_df)

        # Create artificial success_bools to be used by remove_sample_outliers
        success_bools = np.ones(len(dists), dtype=bool)

    return data_df, dists, success_bools, prov_code


def calculate_distances_and_optimize(data_df, offset_bounds):
    """For each sample, perform optimized load balancing.

    This means finding an offset for each sample that minimizes
    the distance metric. The distance metric is the sum of the
    distances from each probe measurement to its median.

    Args:
        data_df: pandas dataframe of floats
        offset_bounds: tuple of floats

    Returns:
        out_df: pandas dataframe with offsets applied
        optim_distances: numpy array with length = num_samples
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
        sample_vals = data_df.loc[:, sample_ind]

        # Calculate optimized distances
        optimization_result = minimize_scalar(distance_function, args=(sample_vals, probe_medians),
                                              method="bounded", bounds=offset_bounds)

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

    # Return the df with offsets applied, the optimized distances,
    # and Boolean array of samples that converged
    return df_with_offsets, optimized_distances, success_bools


def calculate_distances(data_df):
    """For each sample, calculate distance metric.

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
        sample_vals = data_df.loc[:, sample_ind]

        # Calculate unoptimized distances
        dists[sample_ind] = distance_function(0, sample_vals, probe_medians)

    return dists


def distance_function(offset, values, medians):
    """This function calculates the distance metric.

    Args:
        offset: float
        values: numpy array of floats
        medians: numpy array of floats
    Returns:
        dist: float
    """
    dist = sum(np.square((offset + values) - medians))
    return dist


def remove_sample_outliers(data_df, distances, success_bools, sd_sample_outlier_cutoff):
    """Calculate sample outlier cutoff and remove outlier samples.

    Samples are considered outliers if their distance is above the cutoff
    OR if the optimization process didn't converge for them.

    Args:
        data_df: pandas dataframe of floats
        distances: numpy array of floats with length = num_samples
        success_bools: numpy array of bools with length = num_samples
        sd_sample_outlier_cutoff: float
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
                                              sd_sample_outlier_cutoff)

    # Remove samples whose distance metric is greater than the cutoff OR
    # those that didn't converge during optimization
    out_df = data_df.iloc[:, np.logical_and(distances < cutoff, success_bools)]
    assert not out_df.empty, "All samples were filtered out. Try reducing the SD cutoff."
    return out_df

def row_median_normalize(data_df):
    """Subtract the row median from values.

    Args:
        data_df: dataframe of floats
    Returns:
        out_df: dataframe of floats
    """
    out_df = data_df.subtract(data_df.median(axis=1), axis=0)
    return out_df


if __name__ == '__main__':
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args:  {}".format(args))

    main(args)
