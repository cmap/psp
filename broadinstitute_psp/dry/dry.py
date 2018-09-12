"""
dry.py

Performs filtering and normalization of P100 and GCP data.

Converts level 2 to level 3 data. Required input (--in_gct_path) is a path to a
gct file. Output is writing a processed gct file and a pw (plate-well) file
with QC information.

"""

import argparse
import logging
import numpy as np
import os
import pandas as pd
import sys

import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.write_gct as wg
import broadinstitute_psp.utils.qc_gct2pw as gct2pw
import broadinstitute_psp.utils.separate_gct as sg
import broadinstitute_psp.utils.psp_utils as psp_utils
import broadinstitute_psp.utils.setup_logger as setup_logger

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Default output suffixes
DEFAULT_GCT_SUFFIX = ".dry.processed.gct"
DEFAULT_PW_SUFFIX = ".dry.processed.pw"


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required arg
    parser.add_argument("--in_gct_path", "-i", required=True,
                        help="filepath to input gct")

    # Optional args
    parser.add_argument("--out_dir", "-o", default=".",
                        help="path to save directory")
    parser.add_argument("--out_base_name", "-ob", default=None,
                        help="base name of output pw and gct files (default is <INPUT_GCT>.dry.processed")
    parser.add_argument("--psp_config_path", "-p", default="~/psp_production.cfg",
                        help="filepath to PSP config file")
    parser.add_argument("--force_assay", "-f",
                        choices=["GCP", "P100"], default=None,
                        help=("directly specify assay type here " +
                              "(overrides first entry of provenance code)"))
    parser.add_argument("--no_optim", "-no", action="store_true", default=False,
                        help="whether to perform P100 load balancing optimization")

    # Parameters
    parser.add_argument("--sample_frac_cutoff", "-sfc", type=float, default=None,
                        help=("fraction of probes that must be measured in " +
                              "a sample to avoid the sample being filtered " +
                              "out; if None, assay-specific default value " +
                              "from config file will be used"))
    parser.add_argument("--probe_frac_cutoff", "-pfc", type=float, default=None,
                        help=("fraction of samples in which a probe must be " +
                              "measured to avoid the probe being filtered " +
                              "out; if None, assay-specific default value from " +
                              "config file will be used"))
    parser.add_argument("--probe_sd_cutoff", "-psc", type=float, default=None,
                        help=("the SD of a probe must be less than this cutoff " +
                              "to avoid being filtered out; if None, " +
                              "assay-specific default value from " +
                              "config file will be used"))
    parser.add_argument("--dist_sd_cutoff", "-dsc", type=float, default=5,
                        help=("maximum SD for a sample's distance metric " +
                              "before being filtered out"))
    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                        help="increase the number of messages reported")

    return parser


def main(args):
    """THE MAIN METHOD. Filter and normalize data and save result as a gct file.

    Args:
        args (argparse.Namespace object): fields as defined in build_parser()

    Returns:
        out_gct (GCToo object): output gct object
    """
    ### READ GCT AND CONFIG FILE
    (in_gct, assay_type, prov_code, config_io, config_metadata, config_parameters) = (
        read_dry_gct_and_config_file(
            args.in_gct_path, args.psp_config_path, args.force_assay))

    ### LOG TRANSFORM
    (l2x_gct, prov_code) = log_transform_if_needed(
        in_gct, prov_code, config_metadata["log_transform_prov_code_entry"])

    ### INITIAL FILTERING
    (filt_gct, prov_code, post_sample_nan_remaining) = initial_filtering(
        l2x_gct, assay_type, args.sample_frac_cutoff, args.probe_frac_cutoff,
        args.probe_sd_cutoff, config_parameters,
        config_metadata["manual_rejection_field"],
        prov_code, config_metadata["sample_filter_prov_code_entry"],
        config_metadata["manual_probe_reject_prov_code_entry"],
        config_metadata["probe_filter_prov_code_entry"])

    ### HISTONE NORMALIZE (if GCP)
    (hist_norm_gct, prov_code) = gcp_histone_normalize_if_needed(
        filt_gct, assay_type, config_metadata["gcp_normalization_peptide_field"],
        config_metadata["gcp_normalization_peptide_id"], prov_code,
        config_metadata["gcp_histone_prov_code_entry"])

    ### APPLY OFFSETS IF NEEDED (if P100)
    (offset_gct, dists, offsets, prov_code) = p100_calculate_dists_and_apply_offsets_if_needed(
            hist_norm_gct, assay_type, args.no_optim,
            eval(config_parameters["offset_bounds"]), prov_code,
            config_metadata["optimization_prov_code_entry"])

    ### FILTER SAMPLES BY DISTANCE (if P100)
    (filt_dist_gct, out_offsets, post_sample_dist_remaining, prov_code) = p100_filter_samples_by_dist(
            offset_gct, assay_type, offsets, dists,
            args.dist_sd_cutoff, prov_code,
            config_metadata["outlier_sample_filter_prov_code_entry"])

    ### INSERT OFFSETS AND UPDATE PROVENANCE CODE
    out_gct = insert_offsets_and_prov_code(
        filt_dist_gct, out_offsets, config_metadata["offsets_field"], prov_code,
        config_metadata["prov_code_field"], config_metadata["prov_code_delimiter"])

    ### CONFIGURE OUT NAMES
    (out_gct_name, out_pw_name) = configure_out_names(
        args.in_gct_path, args.out_base_name)

    ### WRITE PW FILE OF SAMPLES FILTERED
    write_output_pw(in_gct, post_sample_nan_remaining,
                    post_sample_dist_remaining, offsets,
                    args.out_dir, out_pw_name)

    ### WRITE OUTPUT GCT
    write_output_gct(out_gct, args.out_dir, out_gct_name,
                     config_io["data_null"], config_io["filler_null"])

    return out_gct


# tested #
def read_dry_gct_and_config_file(in_gct_path, config_path, forced_assay_type):
    """ Read gct and config file.

    Uses the utility function read_gct_and_config_file from psp_utils.
    Provenance code is extracted from the col metadata. It must be non-empty
    and the same for all samples. If forced_assay_type is not None,
    assay_type is set to forced_assay_type.

    Args:
        in_gct_path (string): filepath to gct file
        config_path (string): filepath to config file
        forced_assay_type (string, or None)

    Returns:
        gct (GCToo object)
        assay_type (string)
        prov_code (list of strings)
        config_io (dictionary)
        config_metadata (dictionary)
        config_parameters (dictionary)
    """
    # Read gct and config file
    (gct, config_io, config_metadata, config_parameters) = psp_utils.read_gct_and_config_file(in_gct_path, config_path)

    # Extract the plate's provenance code
    prov_code = psp_utils.extract_prov_code(gct.col_metadata_df,
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
    assay_type_out = check_assay_type(assay_type, p100_assay_types, gcp_assay_types)

    return gct, assay_type_out, prov_code, config_io, config_metadata, config_parameters


# tested #
def check_assay_type(assay_type, p100_assays, gcp_assays):
    """Verify that assay is one of the allowed types. Return either P100 or GCP.

    Args:
        assay_type_in (string)
        p100_assays (list of strings)
        gcp_assays (list of strings)

    Returns:
        assay_type_out (string): choices = {"P100", "GCP"}
    """
    if assay_type in p100_assays:
        assay_type_out = "p100"
    elif assay_type in gcp_assays:
        assay_type_out = "gcp"
    else:
        err_msg = ("The assay type is not a recognized P100 or GCP assay. " +
                   "assay_type: {}, p100_assays: {}, gcp_assays: {}")
        logger.error(err_msg.format(assay_type, p100_assays, gcp_assays))
        raise(Exception(err_msg.format(assay_type, p100_assays, gcp_assays)))

    return assay_type_out


# tested #
def configure_out_names(in_gct_path, out_base_name_from_args):
    """If out_base_name_from_args is None, append DEFAULT_GCT_SUFFIX and 
    DEFAULT_PW_SUFFIX to the input gct name to generate the gct and pw file
    respectively.

    Args:
        in_gct_path:
        out_base_name_from_args:

    Returns:
        out_gct_name (file path)
        out_pw_name (file path)

    """    

    if out_base_name_from_args is None:
        input_basename = os.path.basename(in_gct_path)        
        out_gct_name = input_basename + DEFAULT_GCT_SUFFIX
        out_pw_name = input_basename + DEFAULT_PW_SUFFIX
    else:
        input_basename = os.path.basename(out_base_name_from_args)
        out_gct_name = input_basename + DEFAULT_GCT_SUFFIX
        out_pw_name = input_basename + DEFAULT_PW_SUFFIX
        
    return out_gct_name, out_pw_name


# tested #
def log_transform_if_needed(gct, prov_code, prov_code_entry):
    """Perform log2 transformation if it hasn't already been done.

    Args:
        gct (GCToo object)
        prov_code (list of strings)
        prov_code_entry (string)

    Returns:
        out_gct (GCToo object)
        updated_prov_code (list of strings): updated
    """
    # Check if log2 transformation has already occurred
    if prov_code_entry in prov_code:
        logger.info("{} has already occurred.".format(prov_code_entry))
        updated_prov_code = prov_code

        out_gct = gct

    else:
        assert not (gct.data_df < 0).sum().sum(), (
            "data_df should not contain negative values. gct.data_df:\n{}".format(
                gct.data_df))

        out_df = log_transform(gct.data_df, log_base=2)
        updated_prov_code = prov_code + [prov_code_entry]

        # Return new GCToo
        out_gct = GCToo.GCToo(out_df, gct.row_metadata_df, gct.col_metadata_df)

    return out_gct, updated_prov_code


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
def gcp_histone_normalize_if_needed(gct, assay_type, gcp_normalization_peptide_field,
                                    gcp_normalization_peptide_id, prov_code, prov_code_entry):
    """If GCP, normalize to an invariant probe. With the
    addition of H4 probes, we need to use different normalization peptides
    for different rows. This script splits the data_df up according to norm
    peptide and then concats the results together.

    If gcp_normalization_peptide_field is present, extract normalization
    peptides from metadata. If not present, a new column called
    gcp_normalization_peptide_field will be made and populated with
    gcp_normalization_peptide_id.

    Args:
        gct (GCToo object)
        assay_type (string)
        gcp_normalization_peptide_field (string, or None): metadata field in
            GCT that indicates which norm peptide to use for each row
        gcp_normalization_peptide_id (string, or None): if
            gcp_normalization_peptide_field not provided, this peptide_id
            should be used for all rows
        prov_code (list of strings)
        prov_code_entry (string)

    Returns:
        out_gct (GCToo object): data_df and row_metadata_df updated; norm
            peptides have been removed
        updated_prov_code (list of strings)
    """
    if assay_type == "gcp":

        # If needed, create a new column indicating the norm peptide for each probe
        create_norm_peptide_column_if_needed(
            gct.row_metadata_df, gcp_normalization_peptide_field,
            gcp_normalization_peptide_id)

        # Split gct according to norm_peptide
        (split_gcts, list_of_norm_peptides) = sg.separate(gct, gcp_normalization_peptide_field, "row")

        # Make sure that each norm_peptide is assigned to itself
        tmp = gct.row_metadata_df.loc[list_of_norm_peptides, gcp_normalization_peptide_field]
        assert np.array_equal(tmp.index, tmp.values), (
            ("Each normalization peptide must be assigned to itself.\n" +
             "index\t{}\n{}").format(gcp_normalization_peptide_field, tmp))

        # Perform normalization for each subsetted GCT
        list_of_out_dfs = []
        for g, norm_pep in zip(split_gcts, list_of_norm_peptides):
            small_out_df = gcp_histone_normalize(g.data_df, norm_pep)
            list_of_out_dfs.append(small_out_df)

        # Reconcatenate the split up GCTs; note that rows are sorted in output
        out_df = pd.concat(list_of_out_dfs, axis=0).sort_index(axis=0)

        # Remove any cols that have been made all NaN
        out_df = filter_samples_by_nan(out_df, 0)

        # Update gct object and provenance code
        (out_gct, updated_prov_code) = update_metadata_and_prov_code(
            out_df, gct.row_metadata_df, gct.col_metadata_df, prov_code_entry, prov_code)

    else:
        out_gct = gct
        updated_prov_code = prov_code

    return out_gct, updated_prov_code


def create_norm_peptide_column_if_needed(row_meta_df, gcp_normalization_peptide_field, gcp_normalization_peptide_id):
    """If it doesn't already exist, create a new column indicating the norm
    peptide for each probe.

    Modifies row_meta_df in-place.

    Args:
        row_meta_df (pandas df):
        gcp_normalization_peptide_field (string): metadata field in
            GCT that indicates which norm peptide to use for each probe
        gcp_normalization_peptide_id (string): if
            gcp_normalization_peptide_field not in the GCT, a new column will
            be created and filled with just this peptide_id

    Returns:
        None

    """
    # If gcp_normalization_peptide_field is already in the row metadata, we'll
    # just use it
    if gcp_normalization_peptide_field not in row_meta_df.columns.values:

        assert gcp_normalization_peptide_id is not None, (
            "gcp_normalization_peptide_field not present in metadata headers, " +
            "so gcp_normalization_peptide_id must not be None. " +
            "gcp_normalization_peptide_field: {}, " +
            "gcp_normalization_peptide_id: {}").format(
            gcp_normalization_peptide_field, gcp_normalization_peptide_id)

        # Otherwise, create a new column and fill it with gcp_normalization_peptide_id
        row_meta_df[gcp_normalization_peptide_field] = gcp_normalization_peptide_id

    return None


# tested #
def gcp_histone_normalize(data_df, gcp_normalization_peptide_id):
    """Subtract values of gcp_normalization_peptide_id from all the other probes.
    Remove the row of data corresponding to the normalization histone.

    Assumes that all probes should be normalized to the same peptide id.

    Args:
        data_df (pandas df)
        gcp_normalization_peptide_id (string): id

    Returns:
        out_df (pandas df): one or more rows removed
    """
    # Verify that norm peptide is in the data
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
def initial_filtering(gct, assay_type, sample_frac_cutoff, probe_frac_cutoff, probe_sd_cutoff,
                      config_parameters, manual_rejection_field, prov_code,
                      sample_filt_prov_code_entry, manual_reject_prov_code_entry, probe_filt_prov_code_entry):
    """Performs three types of filtering. filter_samples_by_nan removes
    samples that have too many nan values. manual_probe_rejection removes
    probes that were manually labeled for removal. filter_probes_by_nan_and_sd
    removes probes if they have too many nan values OR if the SD for the whole
    row is too high.

    If nan_thresh is None, will use assay-specific default from config file.

    Also return a list of samples remaining after sample filtration.

    Args:
        gct (GCToo object)
        assay_type (string)
        sample_frac_cutoff (float b/w 0 and 1)
        probe_frac_cutoff (float b/w 0 and 1)
        probe_sd_cutoff (float)
        config_parameters (dictionary)
        manual_rejection_field (string)
        prov_code (list of strings)
        sample_filt_prov_code_entry (string)
        manual_reject_prov_code_entry (string)
        probe_filt_prov_code_entry (string)

    Returns:
        out_gct (GCToo object): updated
        post_sample_nan_remaining (list of strings)
        prov_code (list of strings): updated
    """
    # Record what probes we begin with
    initial_probes = list(gct.data_df.index.values)

    # Check assay-specific thresholds
    [sample_frac_cutoff, probe_frac_cutoff, probe_sd_cutoff] = check_assay_specific_thresh(
        assay_type, sample_frac_cutoff, probe_frac_cutoff, probe_sd_cutoff, config_parameters)

    ### FILTER SAMPLES BY NAN
    sample_nan_data_df = filter_samples_by_nan(gct.data_df, sample_frac_cutoff)
    thresh_digit = ("{:.1f}".format(sample_frac_cutoff)).split(".")[1]
    sample_filt_prov_code_entry_formatted = "{}{}".format(sample_filt_prov_code_entry, thresh_digit)

    (out_gct, prov_code) = update_metadata_and_prov_code(
        sample_nan_data_df, gct.row_metadata_df, gct.col_metadata_df, sample_filt_prov_code_entry_formatted, prov_code)

    # Record what samples remain
    post_sample_nan_remaining = list(out_gct.data_df.columns.values)

    ### FILTER MANUALLY REJECTED PROBES
    probe_manual_data_df = manual_probe_rejection(out_gct.data_df, out_gct.row_metadata_df, manual_rejection_field)

    # Check if any probes were actually rejected
    post_probe_manual_remaining_probes = list(probe_manual_data_df.index.values)
    probes_removed = (set(initial_probes) != set(post_probe_manual_remaining_probes))

    # Only update prov code if probes were actually rejected
    if probes_removed:
        (out_gct, prov_code) = update_metadata_and_prov_code(
            probe_manual_data_df, out_gct.row_metadata_df, out_gct.col_metadata_df,
            manual_reject_prov_code_entry, prov_code)

    ### FILTER PROBES BY NAN AND SD
    probe_nan_sd_data_df = filter_probes_by_nan_and_sd(out_gct.data_df, probe_frac_cutoff, probe_sd_cutoff)
    thresh_digit = ("{:.1f}".format(probe_frac_cutoff)).split(".")[1]
    probe_filt_prov_code_entry_formatted = "{}{}".format(probe_filt_prov_code_entry, thresh_digit)
    (out_gct, prov_code) = update_metadata_and_prov_code(
        probe_nan_sd_data_df, out_gct.row_metadata_df, out_gct.col_metadata_df, probe_filt_prov_code_entry_formatted, prov_code)

    return out_gct, prov_code, post_sample_nan_remaining

# tested #
def check_assay_specific_thresh(assay_type, sample_frac_cutoff, probe_frac_cutoff,
                                probe_sd_cutoff, config_parameters):
    """Checks if sample_frac_cutoff, probe_frac_cutoff, or probe_sd_cutoff were
     provided. If not, uses values from config file.

    Args:
        sample_frac_cutoff (float)
        probe_frac_cutoff (float)
        probe_sd_cutoff (float)
        config_parameters (dictionary): contains assay-specific threshold
            values to use if provided thresholds are None

    Returns:
        sample_frac_cutoff_out (float)
        probe_frac_cutoff_out (float)
        probe_sd_cutoff_out (float)

    """
    if sample_frac_cutoff is None:
        if assay_type == "p100":
            sample_frac_cutoff_out = float(config_parameters["p100_sample_frac_cutoff"])
        elif assay_type == "gcp":
            sample_frac_cutoff_out = float(config_parameters["gcp_sample_frac_cutoff"])
    else:
        sample_frac_cutoff_out = sample_frac_cutoff

    if probe_frac_cutoff is None:
        if assay_type == "p100":
            probe_frac_cutoff_out = float(config_parameters["p100_probe_frac_cutoff"])
        elif assay_type == "gcp":
            probe_frac_cutoff_out = float(config_parameters["gcp_probe_frac_cutoff"])
    else:
        probe_frac_cutoff_out = probe_frac_cutoff

    if probe_sd_cutoff is None:
        if assay_type == "p100":
            probe_sd_cutoff_out = float(config_parameters["p100_probe_sd_cutoff"])
        elif assay_type == "gcp":
            probe_sd_cutoff_out = float(config_parameters["gcp_probe_sd_cutoff"])
    else:
        probe_sd_cutoff_out = probe_sd_cutoff

    return sample_frac_cutoff_out, probe_frac_cutoff_out, probe_sd_cutoff_out

# tested #
def filter_samples_by_nan(data_df, sample_frac_cutoff):
    """ Filter out samples whose fraction of probes measured is less than sample_frac_cutoff.

    Args:
        data_df (pandas df)
        sample_frac_cutoff (float b/w 0 and 1)
    Returns:
        out_df (pandas df): potentially smaller than input
    """
    # Number of NaNs per sample
    num_nans = data_df.isnull().sum()

    # Number of rows
    num_rows = data_df.shape[0]

    # Fraction non-NaN per sample
    frac_non_nan_per_sample = 1 - num_nans/num_rows

    # Only return samples with more non-NaN data than sample_frac_cutoff
    out_df = data_df.loc[:, frac_non_nan_per_sample > sample_frac_cutoff]
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
def filter_probes_by_nan_and_sd(data_df, probe_frac_cutoff, probe_sd_cutoff):
    """ Filter out probes whose fraction of samples measured is less than probe_frac_cutoff.
    Also remove probes with standard deviation higher than probe_sd_cutoff.

    Args:
        data_df (pandas df)
        probe_frac_cutoff (float b/w 0 and 1)
        probe_sd_cutoff (float)
    Returns:
        out_df (pandas df): potentially smaller than original df
    """
    # Number of NaNs per probe
    num_nans = data_df.isnull().sum(axis=1)

    # Number of samples
    num_samples = data_df.shape[1]

    # Fraction non-NaN per probe
    frac_non_nans_per_probe = 1 - num_nans/num_samples

    # Probe standard deviations
    probe_sds = data_df.std(axis=1)

    # Only return probes with more non-NaN data than probe_frac_cutoff
    # and lower sd than the cutoff
    probes_to_keep = ((frac_non_nans_per_probe > probe_frac_cutoff) &
                      (probe_sds < probe_sd_cutoff))
    out_df = data_df.loc[probes_to_keep, :]
    assert not out_df.empty, (
        "All probes were filtered out. Try reducing the NaN threshold and/or SD cutoff.")
    return out_df

# tested #
def p100_calculate_dists_and_apply_offsets_if_needed(gct, assay_type, no_optim_bool,
                                                     offset_bounds, prov_code, prov_code_entry):
    """If P100, calculate a distance metric for each sample in order to use it
    later for filtering. The distance metric is how far away each probe is from
    its median value.

    Optimization happens by default for P100, but can be turned off by setting
    no_optim_bool to False. During optimization, an offset is
    calculated for each sample that seeks to minimize the distance metric for
    that sample. The distance is then recalculated for the data after offsets
    have been applied. The provenance code is only updated if optimization occurs.

    Args:
        gct (GCToo object)
        assay_type (string)
        no_optim_bool (bool): if true, optimization will NOT occur
        offset_bounds (tuple of floats): specifies bounds for the offset;
            if outside of the bounds, a warning will be thrown
        prov_code (list of strings)
        prov_code_entry (string)

    Returns:
        out_gct (GCToo object): updated
        dists (numpy array of floats): distance metric for each sample
        offsets (numpy array of floats): offset that was applied to each sample
        prov_code (list of strings): updated

    """
    # P100
    if assay_type == "p100":
        if not no_optim_bool:
            # Perform optimization and return offsets and distances
            (out_df, offsets, dists) = (
                calculate_distances_and_optimize(gct.data_df, offset_bounds))
            (out_gct, prov_code) = update_metadata_and_prov_code(
                out_df, gct.row_metadata_df, gct.col_metadata_df, prov_code_entry, prov_code)
        else:
            # Simply calculate distance metric for each sample
            dists = calculate_distances(gct.data_df)

            # Set offsets to None
            offsets = None
            out_gct = gct

    # GCP
    # N.B. distances are not calculated because filtration by distance doesn't occur
    else:
        dists = None
        offsets = None
        out_gct = gct

    return out_gct, dists, offsets, prov_code


# tested #
def calculate_distances_and_optimize(data_df, offset_bounds):
    """For each sample, perform optimization.

    This means finding an offset for each sample that minimizes
    the distance metric. The distance metric is the sum of the
    distances from each probe measurement to its median.

    N.B. Only uses the non-NaN values from each sample.

    Args:
        data_df (pandas df)
        offset_bounds (tuple of floats): specifies bounds for the offset;
            if outside of the bounds, a warning will be thrown

    Returns:
        out_df (pandas dataframe): with offsets applied
        optimized_offsets (numpy array): constant to add to each sample to
            minimize the distance function; length = num_samples
        optimized_dists (numpy array): distance with optimized offset added;
            length = num_samples
    """
    # Compute probe medians
    orig_probe_medians = data_df.median(axis=1).values

    # Compute offsets for all samples
    optimized_offsets = calculate_offsets_analytically(data_df)

    # Apply the optimized offsets
    df_with_offsets = data_df.add(optimized_offsets, axis=1)

    # Compute distances after applying offsets
    optimized_dists = [distance_function(df_with_offsets.loc[:, col].values, orig_probe_medians) for col in df_with_offsets]

    # Report which samples had offsets greater than the bounds
    assert offset_bounds[0] < offset_bounds[1]
    offsets_outside_of_bounds_bool_array = (optimized_offsets < offset_bounds[0]) | (optimized_offsets > offset_bounds[1])
    if offsets_outside_of_bounds_bool_array.any():
        offsets_outside_of_bounds = df_with_offsets.columns[offsets_outside_of_bounds_bool_array].values
        offsets_outside_of_bounds_vals = optimized_offsets[offsets_outside_of_bounds_bool_array].values
        logger.warning((
            "The following samples have optimized offsets outside the " +
            "offset_bounds.\n{}").format(
            zip(offsets_outside_of_bounds, offsets_outside_of_bounds_vals)))

    # Return the df with offsets applied, the offsets, and the optimized dists
    return df_with_offsets, optimized_offsets, optimized_dists


def calculate_offsets_analytically(data_df):
    """ Calculate offsets analytically.

    Equation for calculating offsets:
    o = (P - s) / N

    o is the vector of offsets
    P is the scalar value of the sum of probe medians
    s is the vector of sums of each sample's values
    N is the scalar of the number of probes

    One subtlety is that if there is missing data for a sample, the
    corresponding probe median must also be omitted.

    Args:
        data_df (pandas df)

    Returns:
        optimized_offsets (pandas series)

    """
    # Calculate probe medians
    probe_medians = data_df.median(axis=1)

    # Initialize output
    optimized_offsets = pd.Series(np.full((data_df.shape[1]), np.nan), index=data_df.columns)

    for col in data_df:

        # Determine which probes have missing data
        bool_array_of_not_nan = pd.notnull(data_df.loc[:, col])

        # Take the sum of only those probe medians
        this_sum_of_probe_medians = probe_medians.loc[bool_array_of_not_nan].sum()

        # Calculate the sum of values for this sample (NaN excluded by default)
        this_sum_of_sample_values = data_df.loc[:, col].sum()

        # Count number of non-NaN probes
        num_probes = bool_array_of_not_nan.sum()

        # Calculate offset
        this_offset = (this_sum_of_probe_medians - this_sum_of_sample_values) / num_probes

        # Append offset
        optimized_offsets[col] = this_offset

    return optimized_offsets


# tested #
def calculate_distances(data_df):
    """For each sample, calculate distance metric.
    N.B. Only uses the non-NaN values from each sample.

    Args:
        data_df (pandas df)
    Returns:
        dists (numpy array of floats): length = num_samples
    """
    # Determine the median value for each probe
    probe_medians = data_df.median(axis=1)

    # Calculate distances
    dists = [distance_function(data_df.loc[:, col].values, probe_medians) for col in data_df]

    return dists

# tested #
def distance_function(values, medians):
    """This function calculates the distance metric.
    N.B. Only uses the non-NaN values.

    dist = sum( (s - m)^2 )

    s is the vector of sample values
    m is the vector of probe medians

    Args:
        values (numpy array of floats)
        medians (numpy array of floats)
    Returns:
        dist (float)
    """
    non_nan_idx = ~np.isnan(values)
    assert np.size(non_nan_idx) != 0, "All values in this sample are NaN!"

    non_nan_values = values[non_nan_idx]
    non_nan_medians = medians[non_nan_idx]
    dist = sum(np.square(non_nan_values - non_nan_medians))
    return dist

# tested #
def p100_filter_samples_by_dist(gct, assay_type, offsets, dists,
                                dist_sd_cutoff, prov_code, prov_code_entry):
    """If P100, filter out samples whose distance metric is above some threshold.
    Also remove samples for which optimization did not converge.

    Also return a list of samples that remain after filtration.

    N.B. Offsets are only passed to this function so that offsets of filtered
    out samples can be removed.

    Args:
        gct (GCToo object)
        assay_type (string)
        offsets (pandas series of floats, or None): offset that was applied to each sample
        dists (numpy array of floats, or None): distance metric for each sample
        dist_sd_cutoff (float): maximum SD for a sample's distance metric before being filtered out
        prov_code (list of strings)
        prov_code_entry (string)

    Returns:
        gct (GCToo object): updated
        out_offsets (numpy array of floats, or None): only returned for samples that were not removed
        post_sample_dist_remaining (list of strings): if GCP, then None
        prov_code (list of strings): updated

    """
    # P100
    if assay_type == "p100":

        (out_df, out_offsets) = remove_sample_outliers(gct.data_df, offsets, dists, dist_sd_cutoff)
        prov_code_entry_formatted = "{}{:.0f}".format(prov_code_entry, dist_sd_cutoff)
        (gct, prov_code) = update_metadata_and_prov_code(
            out_df, gct.row_metadata_df, gct.col_metadata_df, prov_code_entry_formatted, prov_code)

         # Record what samples remain
        post_sample_dist_remaining = list(gct.data_df.columns.values)

    # GCP
    else:
        out_offsets = None
        post_sample_dist_remaining = None

    return gct, out_offsets, post_sample_dist_remaining, prov_code


# tested #
def remove_sample_outliers(data_df, offsets, distances, dist_sd_cutoff):
    """Calculate distance cutoff for outlier samples and remove outliers.

    Samples are considered outliers if their distance is above the cutoff.
    Return only the offsets of samples that were not removed.

    Args:
        data_df (pandas df)
        offsets (pandas series of floats): length = num_samples
        distances (numpy array of floats): length = num_samples
        dist_sd_cutoff (float): maximum SD for a sample's distance metric before being filtered out

    Returns:
        out_df (pandas df): potentially smaller than original df
        out_offsets (numpy array of floats): length = num_remaining_samples
    """
    assert len(distances) == data_df.shape[1], (
        "len(distances): {} does not equal data_df.shape[1]: {}").format(len(distances), data_df.shape[1])

    # Calculate the distance cutoff
    cutoff = np.mean(distances) + np.multiply(np.std(distances), dist_sd_cutoff)

    # Determine what samples to keep
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
        gct.col_metadata_df[offsets_field] = offsets

    # Convert provenance code to delimiter separated string
    prov_code_str = prov_code_delimiter.join(prov_code)

    # Update the provenance code in col_metadata_df
    gct.col_metadata_df.loc[:, prov_code_field] = prov_code_str

    return gct

# tested #
def write_output_pw(gct, post_sample_nan_remaining, post_sample_dist_remaining,
                    offsets, out_dir, out_pw_name):
    """Save to a .pw file a record of what samples were filtered out and all
    optimization offsets that were computed.

    Args:
        gct (GCToo object): original input gct
        post_sample_nan_remaining (list of strings):
        post_sample_dist_remaining (list of strings): if GCP, this will be None
        offsets (pandas series of floats): if GCP, this will be None
        out_dir (filepath): where to save output file
        out_pw_name (string): what to call output file

    Returns:
        .pw file called out_pw_name

    """
    # Extract plate and well names
    [plate_names, well_names] = gct2pw.extract_plate_and_well_names(
        gct.col_metadata_df, "det_plate", "det_well")

    # Create boolean lists of length = number of original samples indicating if
    # sample was filtered at a certain step
    original_samples = list(gct.data_df.columns.values)
    post_sample_nan_bools = [sample in post_sample_nan_remaining for sample in original_samples]

    if post_sample_dist_remaining is not None:
        post_sample_dist_bools = [sample in post_sample_dist_remaining for sample in original_samples]

        if offsets is not None:
            offsets_to_insert = [offsets.loc[sample] if sample in offsets.index else np.nan for sample in original_samples]
            # Assemble plate_names, well_names, bool arrays, and offsets into a df
            out_df = gct2pw.assemble_output_df(
                plate_names, well_names, {
                    "remains_after_poor_coverage_filtration": post_sample_nan_bools,
                    "remains_after_outlier_removal":post_sample_dist_bools,
                    "optimization_offset":offsets_to_insert})
        else:
            # Assemble plate_names, well_names, and bool arrays into a df
            out_df = gct2pw.assemble_output_df(
                plate_names, well_names, {
                    "remains_after_poor_coverage_filtration": post_sample_nan_bools,
                    "remains_after_outlier_removal":post_sample_dist_bools})
    else:
        # Assemble plate_names, well_names, and single bool array into a df
        out_df = gct2pw.assemble_output_df(
            plate_names, well_names,{
                "remains_after_poor_coverage_filtration":post_sample_nan_bools})

    # Write to pw file
    full_out_pw_name = os.path.join(out_dir, out_pw_name)
    out_df.to_csv(full_out_pw_name, sep="\t", na_rep="NaN", index=False)
    logger.info("PW file written to {}".format(full_out_pw_name))


def write_output_gct(gct, out_dir, out_gct_name, data_null, filler_null):
    """Write output gct file.

    Args:
        gct (GCToo object)
        out_dir (string): path to save directory
        out_gct_name (string): name of output gct
        data_null (string): string with which to represent NaN in data
        filler_null (string): string with which to fill the empty top-left quadrant in the output gct

    Returns:
        None

    """
    out_fname = os.path.join(out_dir, out_gct_name)
    wg.write(gct, out_fname, data_null=data_null, filler_null=filler_null, data_float_format=None)


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
# TODO(lev): refactor to use subset_gctoo
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

    main(args)
