import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse
import logging
import sys

import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.write_gct as write_gct
import cmapPy.pandasGEXpress.setup_GCToo_logger as setup_logger
import cmapPy.pandasGEXpress.GCToo as GCToo


logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--in_gct_path", "-i",
                        help="path to GCT to be re-normalized")
    group.add_argument("--in_gct", help="GCT to be re-normalized" )

    # Optional args
    parser.add_argument("--out_name", "-o", type=str, default="renormed",
                        help="what to name the output file " 
                            "(appropriate GCT suffix added automatically)")
    parser.add_argument("--slope_cutoff", "-sc", type=float, default=0.2,
                        help="default slope value to be re-normalized")
    parser.add_argument("--write_gct", "-gct", action="store_true", default=False,
                       help="whether or not to create a GCT of output")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to print a bunch of output")

    return parser



def continuous_renormalization(args):

    # Read in GCT, if path provided, and make deep copies of all DataFrames
    if args.in_gct_path:
        gct = parse.parse(args.in_gct_path)
    else:
        gct = args.in_gct
    data_df = gct.data_df.copy(deep=True)
    row_metadata_df = gct.row_metadata_df.copy(deep=True)
    col_metadata_df = gct.col_metadata_df.copy(deep=True)

    # Remove rows that are all NA values
    data_df = data_df.loc[(data_df.isnull().apply(np.sum, axis=1) < data_df.shape[1]), :]

    # Pull out enrichment scores from column metadata dataframe
    enrichment_scores = col_metadata_df.loc[:, "det_well_enrichment_score"].copy(deep=True)

    # Calculate limit as x approaches 1 for non-median normalized data
    pep_y_offsets = calculate_y_offsets(data_df, enrichment_scores)
    
    # Calculate the fit parameters
    fit_parameters = calculate_fit(data_df, enrichment_scores, pep_y_offsets)

    # Annotate which rows will be renormalized based on slope_cutoff argument (default 0.2)
    row_metadata_df["is_log_renormed"] = is_log_renormed(fit_parameters.loc[:, "deg1"].apply(get_slope),
                                                         args.slope_cutoff)
    
    # Calculate the offset matrix
    offset_mat = calculate_peptide_sample_offsets(data_df, row_metadata_df, enrichment_scores, fit_parameters,
                                                  pep_y_offsets)

    # Calculate the output DataFrame
    out_data_df = calculate_out_matrix(data_df, offset_mat)

    # Add the 'renorm_correction' metadata field with total sample offset values
    col_metadata_df["renorm_correction"] = calculate_total_sample_offsets(offset_mat)

    # Output
    if args.write_gct:
        write_gct.write(GCToo.GCToo(data_df=out_data_df,
                                col_metadata_df=col_metadata_df,
                                row_metadata_df=row_metadata_df),
                        args.out_name)
    else:
        return GCToo.GCToo(data_df=out_data_df,
                           col_metadata_df=col_metadata_df,
                           row_metadata_df=row_metadata_df)


def calculate_y_offsets(data_df, enrichment_scores, top_fraction=0.2):
    """ Calculate the y-axis offsets by taking the median of the
    top fraction (default 20%) of enrichment scores per peptide

    Args :
        data_df (DataFrame)
        enrichment_scores (DataFrame)
        top_fraction (float) - fraction of enrichment scores to use
            for median calculation

    Return :
        pep_y_offsets (Series) - offsets with index of data_df

    """

    # Calculate number of samples needs for top fraction and set up output
    n_samples_for_median = int(data_df.shape[1] * top_fraction)
    median_enrichment_score_offsets = pd.Series(index=data_df.index)

    for row in median_enrichment_score_offsets.index:

        # Filter out NA samples
        na_samples = data_df.loc[row, :].notnull()
        na_filtered_peptides = data_df.loc[row, na_samples]
        na_filtered_enrichment_scores = enrichment_scores.loc[na_samples]

        # Pull out N largest enrichment scores from filtered scores
        highest_es_samples = na_filtered_enrichment_scores.nlargest(n_samples_for_median).index

        # Pull out dataframe of filtered peptides and calculate median
        # for N largest enrichment scores, add to output
        df_slice_for_median = na_filtered_peptides.loc[highest_es_samples]
        median_enrichment_score_offsets.loc[row] = df_slice_for_median.median()
    
    return median_enrichment_score_offsets


def calculate_fit(data_df, es, pep_y_offsets):

    fit_params = pd.DataFrame(index=data_df.index,
                              columns=["deg1", "log"])

    for peptide_name, _ in fit_params.iterrows():

        for degree, _ in fit_params.iteritems():
            (fit_params.loc[peptide_name, degree], y) = fit_single_peptide(data_df, es, peptide_name,
                                                                         degree, pep_y_offsets)

    return fit_params


def fit_single_peptide(data_df, enrichment_scores, peptide_name, degree, peptide_y_offsets):
    """
    Args:
        data_df (DataFrame)
        enrichment_scores (DataFrame)
        peptide_name (string) - name of peptide used pull out data and y offsets
        degree (int) - specifies fit calculation
        peptide_y_offsets (Series) - y-axis offsets based on median top fraction enrichment scores

    Returns:
        fit
        y (DataFrame) - slice of DF containing y values
    """

    y = data_df.loc[peptide_name, :]
    y_offset = peptide_y_offsets.loc[peptide_name]
    samples_used = y.notnull()

    if degree == "log":
        fit = fit_log(enrichment_scores[samples_used], y[samples_used], y_offset)
    else:
        fit = np.polyfit(enrichment_scores[samples_used],
                         y[samples_used],
                         deg=1)
    return (fit, y)


def fit_log(x, quant_values, y_offset):

    yint_min = np.min(quant_values)
    yint_max = np.max(quant_values)
    quant_values = [np.mean(s) for s in quant_values]

    logist = calculate_logistic_function(y_offset)

    (a, b), covariance = curve_fit(logist,
                                   x,
                                   quant_values,
                                   bounds=((2 * yint_min, 0),
                                           (2 * yint_max, np.inf)))
    return (a, b)


def is_log_renormed(slope_vector, slope_cutoff):
    abs_slope_vector = abs(slope_vector)
    return abs_slope_vector > slope_cutoff


def get_slope(fit):
    return fit[0]


def calculate_peptide_sample_offsets(data_df, row_metadata_df, enrichment_scores,
                                     fit_parameters, peptide_y_offsets):
    """
    Args:
        data_df (DataFrame) - data
        row_metadata_df (DataFrame) - row metadata
        enrichment_scores (DataFrame) - enrichment scores
        fit_parameters
        peptide_y_offsets (Series) - y-axis offsets based on median top fraction enrichment scores
    Returns :
        offset_matrix (DataFrame)
    """
    to_be_renormed = data_df[row_metadata_df.is_log_renormed].index
    offset_matrix = pd.DataFrame(0,
                              index=data_df.index,
                              columns=data_df.columns)
    for peptide in to_be_renormed:
        y_offset = peptide_y_offsets.loc[peptide]
        offset_matrix.loc[peptide, :] = make_y(enrichment_scores,
                                               fit_parameters.loc[peptide, "log"],
                                               y_offset)
    return offset_matrix


def make_y(x, model, y_offset=None):
    # Model is linear or flat
    if y_offset is None:
        y = np.zeros(x.shape)
        degree = 1
        for i in xrange(degree):
            y += (x ** (degree - i - 1)) * model[i]
    # Model is logistic
    else:
        logistic_model = calculate_logistic_function(y_offset)
        y = logistic_model(x, *model)
    return y


def calculate_logistic_function(y_offset):
    """ Create logistic function that integrates peptide-specific y-axis offsets """

    def logist(x, a, b):
        return (a / (1 + np.exp(b * x))) + y_offset

    return logist

def calculate_out_matrix(df, offset_mat):
    normed_df = df.sub(offset_mat)
    return normed_df


def calculate_total_sample_offsets(offset_mat):
    abs_total_offset = abs(offset_mat).apply(np.sum)
    return abs_total_offset



def plot_pep(data_df, es, fit_params, pep_name, pep_y_offset, degree=None):
    plt.figure()
    plt.scatter(es,
                data_df.loc[pep_name, :])
    plt.ylabel("Quant Value")
    plt.xlabel("Enrichment Score")

    if degree is not None:
        x = np.linspace(np.min(es), 1, 101)
        plt.plot(x,
                 make_y(x,
                        fit_params.loc[pep_name, degree],
                        pep_y_offset))
    plt.title(pep_name)
    return plt



if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args: {}".format(args))

    continuous_renormalization(args)
