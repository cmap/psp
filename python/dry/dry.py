import logging
import utils.setup_logger as setup_logger
import ConfigParser
import argparse
import sys

import numpy as np
from scipy.optimize import minimize_scalar
import pandas as pd
import in_out.parse_gctoo as parse_gct


logger = logging.getLogger(setup_logger.LOGGER_NAME)

# TO-DO(lev): update functions to retain intermediate results?
# (final gct doesn't have them...)

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
    parser.add_argument("-optim", action="store_true",
                        help="whether to perform load balancing optimization")
    parser.add_argument("-log2", action="store_true",
                        help="whether to perform log2 normalization")
    parser.add_argument("-sample_pct_cutoff", type=float, default=0.3,
                        help=("what percent of sample data can be NaN " +
                              "before being filtered out"))
    parser.add_argument("-probe_pct_cutoff", type=float, default=0.3,
                        help=("what percent of probe data can be NaN " +
                              "before being filtered out"))
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

    return parser


def main(args):
    """The main method. Writes processed gct to a file and produces QC figures.

    Args:
        args: object with fields as defined in the build_parser() method

    Returns:
        out_gct: gct object that has undergone processing
    """
    # Extract what values to consider NaN from config file
    # N.B. eval used to convert string from config file to list
    PSP_nan_values = eval(configParser.get("io", "nan_values"))

    # filepath = "/Users/lev/code/PSP/python/functional_tests/LJP.gct"
    # args.gct_path = "/Users/lev/code/PSP/python/functional_tests/p100_prm_plate29_3H.gct"

    # Parse the gct file and return gct object
    gct = parse_gct.parse(args.gct_path, nan_values=PSP_nan_values)

    # Extract the plate's provenance code and assay type
    prov_code = extract_prov_code(gct.col_metadata_df)
    assay_type = prov_code[0]

    # Make sure assay_type is one of the allowed values
    allowed_assay_types = eval(configParser.get("metadata", "allowed_assay_types"))
    if assay_type not in allowed_assay_types:
        err_msg = ("The assay type is not in the allowed assay types. " +
                   "assay_type: {}, allowed_assay_types: {}")
        logger.error(err_msg.format(assay_type, allowed_assay_types))
        raise(Exception(err_msg.format(assay_type, allowed_assay_types)))

    # Check if log2 transformation has already occurred
    if "L2X" in prov_code:
        logger.info("L2X has already occurred.")
    else:
        # Perform log2 transformation and update provenance code
        gct.data_df = log_transform(gct.data_df)
        prov_code = update_prov_code("L2X", prov_code)

    # TO-DO(lev): If GCP, perform GCP histone normalization

    # Filter samples based on quantity of NaN
    gct.data_df = filter_samples(gct.data_df,
                                 sample_pct_cutoff=args.sample_pct_cutoff)

    # Filter probes that were manually designated for rejection
    (gct.data_df,
     gct.row_metadata_df) = manual_probe_rejection(gct.data_df,
                                                   gct.row_metadata_df)

    # Filter probes based on quantity of NaN and standard deviation
    gct.data_df = filter_probes(gct.data_df,
                                probe_pct_cutoff=args.probe_pct_cutoff,
                                probe_sd_cutoff=args.probe_sd_cutoff)

    # Perform load balancing (only for P100)
    if args.optim and assay_type in ["PR1", "DIA1"]:
        # 1) optimize sample balance
        # 2) apply offsets
        # 3) remove samples that didn't converge and those with distances too great
        (gct.data_df, distances, success_bools) = optimize_sample_balance(df, offset_bounds=(-7, 7))
    else:
        # 1) calculate distances
        # 2) remove samples with distances too great
        (df, distances, success_bools) = calculate_unoptim_distances(df)

    # remove_sample_outliers(df, distances, sd_sample_outlier_cutoff=3)

    # Row median normalize

    # Probe-group specific normalization

    # Save output QC figures

    # Reinsert prov code
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
        return prov_code
    else:
        # TO-DO(lev): should this be different?
        err_msg = ("All columns should have the same provenance code, " +
                   "but actually np.unique(prov_code_list_series.values) = {}")
        logger.error(err_msg.format(np.unique(prov_code_list_series.values)))
        raise(Exception(err_msg.format(np.unique(prov_code_list_series.values))))


def update_prov_code(new_entry, existing_prov_code):
    """Append entry to provenance code.

    Also check that the new entry is an allowed provenance code value.

    Args:
        new_entry: string
        existing_prov_code: list of strings
    Returns:
        updated_prov_code : list of strings
    """
    # TO-DO(lev): check that new_entry is an allowed prov code value

    updated_prov_code = existing_prov_code.extend(new_entry)
    return updated_prov_code


def log_transform(data_df, log_base=2):
    """Take the log of each value in a pandas dataframe.

    Args:
        data_df: pandas dataframe of floats
        log_base: the base of the logarithm (default = 2).

    Returns:
        transformed_data_df: pandas dataframe with log transformed values
    """
    # Data should be a np.ndarray
    assert isinstance(data_df, pd.DataFrame), ("data_df must be a pandas dataframe, " +
                                               "not type(data_df): {}").format(type(data_df))

    # Replace 0 with np.nan
    data_df.replace(0, np.nan, inplace=True)

    # Numpy operations work fine on dataframes
    transformed_data_df = np.log(data_df) / np.log(log_base)
    return transformed_data_df


def manual_probe_rejection(data_df, row_metadata_df):
    """Remove probes that were manually selected for rejection.

    Args:
        data_df: pandas dataframe of floats
        row_metadata_df: pandas dataframe of strings

    Returns:
        out_data_df: pandas dataframe of floats (potentially smaller than input)
        out_meta_df: pandas dataframe of strings (potentially smaller than input)
    """
    # Extract "pr_probe_suitability_manual" metadata field
    keep_probe_str = row_metadata_df.loc[:, "pr_probe_suitability_manual"]

    # Convert strings to booleans
    keep_probe_bool = (keep_probe_str == "TRUE")

    # Check that the list of good probes is not empty
    assert keep_probe_bool.any(), ("No probes were marked TRUE " +
                                   "(i.e. suitable).\nrow_metadata_df.loc[:, " +
                                   "'pr_probe_suitability_manual']: \n{}").format(keep_probe_str)

    # Return the good probes
    out_meta_df = row_metadata_df[keep_probe_bool.values]
    out_data_df = data_df[keep_probe_bool.values]
    return out_data_df, out_meta_df


def filter_samples(data_df, sample_pct_cutoff=0.3):
    """Remove samples (i.e. columns) with more than sample_pct_cutoff NaN values.

    Args:
        data_df: pandas dataframe of floats
        sample_pct_cutoff: float from 0 to 1 (default = 0.3)
    Returns:
        out_df: pandas dataframe (potentially smaller than original df)
    """
    # Input should be a pandas dataframe
    assert isinstance(data_df, pd.DataFrame), ("data_df must be a pandas dataframe, " +
                                               "not type(data_df): {}").format(type(data_df))

    # Number of NaNs per sample
    num_nans = data_df.isnull().sum()

    # Number of rows
    num_rows = data_df.shape[0]

    # Percent NaNs per sample
    pct_nans_per_sample = num_nans / num_rows

    # Only return samples with fewer % of NaNs than the cutoff
    out_df = data_df.loc[:, pct_nans_per_sample < sample_pct_cutoff]
    return out_df


def filter_probes(data_df, probe_pct_cutoff=0.3, probe_sd_cutoff=3):
    """Remove probes (i.e. rows) with more than probe_pct_cutoff NaN values.

    Also remove probes with standard deviation higher than probe_sd_cutoff.

    Args:
        data_df: pandas dataframe of floats
        probe_pct_cutoff: float from 0-1 (default = 0.3)
        probe_sd_cutoff: float (default = 3)
    Returns:
        out_df: pandas dataframe (potentially smaller than original df)
    """
    # Input should be a pandas dataframe
    assert isinstance(data_df, pd.DataFrame), ("data_df must be a pandas dataframe, " +
                                               "not type(data_df): {}").format(type(data_df))

    # Number of NaNs per probe
    num_nans = data_df.isnull().sum(axis=1)

    # Number of samples
    num_samples = data_df.shape[1]

    # Percent NaNs per probe
    pct_nans_per_probe = num_nans / num_samples

    # Probe standard deviations
    probe_sds = data_df.std(axis=1)

    # Only return probes with fewer % of NaNs than the cutoff
    # and lower sd than the cutoff
    probes_to_keep = ((pct_nans_per_probe < probe_pct_cutoff) &
                      (probe_sds < probe_sd_cutoff))
    out_df = data_df.loc[probes_to_keep, :]
    return out_df










### FIGURE OUT THE BEST WAY TO DO THIS.












def calculate_distances(data_df, optim=False, offset_bounds=(-7, 7)):
    """For each sample, calculate distance from probe medians.
    If optim=True, perform optimized load balancing.

    Args:
        data_df: pandas dataframe of floats
        optim: bool for whether to optimize the calculation of sample offset
        offset_bounds: tuple of floats, default = (-7, 7)

    Returns:
        distances: numpy array with length = num_samples
        optimized_offsets:
        success_bools: numpy array of bools indicating if optimization converged
            with length = num_samples (automatically all True if optim=False)
    """
    # Determine the median value for each probe
    probe_medians = data_df.median(axis=1)

    if optim:
        # Initialize optimization outputs
        num_samples = data_df.shape[1]
        optimized_offsets = np.zeros(num_samples, dtype=float)
        optimized_distances = np.zeros(num_samples, dtype=float)
        unoptimized_distances = np.zeros(num_samples, dtype=float)
        success_bools = np.zeros(num_samples, dtype=bool)
        # N.B. 0 is the same as False if you specify that dtype is boolean

    # For each sample, perform optimization
    for sample_ind in range(num_samples):
        sample_vals = data_df.loc[:, sample_ind]

        # Calculate unoptimized distances
        unoptimized_distances[sample_ind] = function_to_optimize(0,
                                                                 sample_vals,
                                                                 probe_medians)

        # Calculate optimized distances
        optimization_result = minimize_scalar(function_to_optimize,
                                              args=(sample_vals, probe_medians),
                                              method="bounded",
                                              bounds=offset_bounds)

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


def function_to_optimize(offset, values, medians):
    """This function calculates a distance metric.

    Args:
        offset: float
        values: numpy array of floats
        medians: numpy array of floats
    Returns:
        out: float
    """
    out = sum(np.square((offset + values) - medians))
    return out


def remove_sample_outliers(data_df, distances, success_bools, sd_sample_outlier_cutoff):
    """Calculate sample outlier cutoff and remove outlier samples.

    Samples are considered outliers if their distance is above the cutoff
    or if the optimization process didn't converge for them.

    Args:
        data_df: pandas dataframe of floats
        distances: numpy array of floats with length = num_samples
        success_bools: numpy array of bools with length = num_samples
        sd_sample_outlier_cutoff: float

    Returns:
        out_df: pandas dataframe (potentially smaller than original df)
    """
    assert len(distances) == data_df.shape[1], ("len(distances): {} does not equal " +
                                                "data_df.shape[1]: {}").format(len(distances), data_df.shape[1])
    assert len(success_bools) == data_df.shape[1], ("len(success_bools): {} does not equal " +
                                                    "data_df.shape[1]: {}").format(len(success_bools), data_df.shape[1])

    # Calculate the distance cutoff
    cutoff = np.mean(distances) + np.multiply(np.std(distances),
                                              sd_sample_outlier_cutoff)

    # Remove samples whose distance metric is greater than the cutoff or
    # those that didn't converge during optimization
    out_df = data_df.iloc[:, np.logical_and(distances < cutoff, success_bools)]
    return out_df


if __name__ == '__main__':
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args:  {}".format(args))

    main(args)
