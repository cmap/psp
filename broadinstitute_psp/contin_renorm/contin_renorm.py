import pandas as pd
import numpy as np
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
    parser.add_argument("--in_gct_path", "-i", required=True,
        help="fields to remove from the common metadata")

    # Optional args
    parser.add_argument("--out_name", "-o", type=str, default="renormed",
        help="what to name the output file " +
             "(appropriate gct suffix added automatically)")
    parser.add_argument("--slope_cutoff", "-sc", default=0.2,
        help="default slope value to be renormalized")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
        help="whether to print a bunch of output")

    return parser



def contin_renorm_main(args):
    """ Separate method from main() in order to make testing easier and to
    enable command-line access. """

    # print in_gct.data_df.isnull().apply(np.sum, axis=1)

    gct = parse.parse(args.in_gct_path)
    data_df = gct.data_df.copy(deep=True)
    row_metadata_df = gct.row_metadata_df.copy(deep=True)
    col_metadata_df = gct.col_metadata_df.copy(deep=True)

    ### hack to remove rows that are all NA values

    data_df = data_df.loc[(data_df.isnull().apply(np.sum, axis=1) < data_df.shape[1]), :]

    # enrichment_score
    es = col_metadata_df.loc[:, "det_well_enrichment_score"].copy(deep=True)

    # calculate lim as x approaches 1 for non-median normalized data
    pep_y_offsets = calc_y_offsets(data_df, es)
    
    # calculate the fit_params
    fit_params = calc_fit(data_df, es, pep_y_offsets)

    # annotate which need to be renormalized
    row_metadata_df["is_log_renormed"] = is_log_renormed(fit_params.loc[:, "deg1"].apply(get_slope),
                                                         args.slope_cutoff)
    
    # calculate the offset matrix
    offset_mat = calc_pep_samp_offsets(data_df, row_metadata_df, es, fit_params,
                                       pep_y_offsets)

    # calculate the output data
    out_data_df = calc_out_mat(data_df, offset_mat)

    # add the metadata field
    col_metadata_df["renorm_correction"] = calc_tot_samp_offsets(offset_mat)

    #return GCToo.GCToo(data_df=out_data_df,
    #                    col_metadata_df=col_metadata_df,
    #                    row_metadata_df=row_metadata_df)

    # write the file
    write_gct.write(GCToo.GCToo(data_df=out_data_df,
                                col_metadata_df=col_metadata_df,
                                row_metadata_df=row_metadata_df),
                    args.out_name)


def calc_y_offsets(data_df, es, top_frac=0.2):
    # make the scatter plots
    # find the median value for the highest enrichment scores
    # return pep_y_offsets a Series with index of data_df that contains
    # all of the offets


    # take the median of the top 20%
    n_samps_for_median = int(data_df.shape[1] * top_frac)
    median_of_highest = pd.Series(index=data_df.index)
    for row in median_of_highest.index:

        na_samps = data_df.loc[row, :].notnull()
        na_filt_pep = data_df.loc[row, na_samps]
        na_filt_es = es.loc[na_samps]

        highest_es_samps = na_filt_es.nlargest(n_samps_for_median).index

        df_slice_for_median = na_filt_pep.loc[highest_es_samps]

        median_of_highest.loc[row] = df_slice_for_median.median()
    
    return median_of_highest


def is_log_renormed(slope_vector, slope_cutoff):
    abs_slope_vector = abs(slope_vector)
    return abs_slope_vector > slope_cutoff


def calc_out_mat(df, offset_mat):
    normed_df = df.sub(offset_mat)
    return normed_df


def calc_tot_samp_offsets(offset_mat):
    abs_tot_offset = abs(offset_mat).apply(np.sum)
    return abs_tot_offset


def calc_pep_samp_offsets(data_df, row_metadata_df, es,
                          fit_params, pep_y_offsets):
    to_be_renormed = data_df[row_metadata_df.is_log_renormed].index
    offset_mat = pd.DataFrame(0,
                              index=data_df.index,
                              columns=data_df.columns)
    for pep in to_be_renormed:
        y_offset = pep_y_offsets.loc[pep]
        offset_mat.loc[pep, :] = make_y(es,
                                        fit_params.loc[pep, "log"],
                                        y_offset)
    return offset_mat


def calc_fit(data_df, es, pep_y_offsets):
    fit_params = pd.DataFrame(index=data_df.index,
                              columns=["deg1", "log"])
    for pep_name, pep_vals in fit_params.iterrows():
        for deg, deg_vals in fit_params.iteritems():
            fit_params.loc[pep_name, deg], y = fit_pep(data_df, es, pep_name,
                                                       deg, pep_y_offsets)
    
    return fit_params


def fit_pep(data_df, es, pep_name, deg, pep_y_offsets):

    y = data_df.loc[pep_name, :]
    y_offset = pep_y_offsets.loc[pep_name]
    samps_used = y.notnull()
    if deg == "log":
        fit = fit_log(es[samps_used], y[samps_used], y_offset)
    else:
        fit = np.polyfit(es[samps_used],
                         y[samps_used],
                         deg=1)
    return fit, y
    

def plot_pep(data_df, es, fit_params, pep_name, pep_y_offset, deg=None):
    plt.figure()
    plt.scatter(es,
                data_df.loc[pep_name, :])
    plt.ylabel("Quant Value")
    plt.xlabel("Enrichment Score")

    if deg is not None:
        x = np.linspace(np.min(es), 1, 101)
        plt.plot(x,
                 make_y(x,
                        fit_params.loc[pep_name, deg],
                        pep_y_offset))
    plt.title(pep_name)
    return plt

    
def make_y(x, model, y_offset=None):
    # model is linear or flat
    if y_offset is None:
        y = np.zeros(x.shape)
        deg = 1
        for i in xrange(deg):
            y += (x ** (deg - i - 1)) * model[i]
    # model is logistic
    else:
        logist_model = calc_logist(y_offset)
        y = logist_model(x, *model)
    return y


def calc_logist(y_offset):
    
    def logist(x, a, b):
        return (a / (1 + np.exp(b * x))) + y_offset

    return logist

#def logist(x, a, b):
#    return a / (1 + np.exp(b * x)) + 2

def fit_log(x, quant_vals, y_offset):

    yint_min = np.min(quant_vals)
    yint_max = np.max(quant_vals)
    quant_vals = [np.mean(s) for s in quant_vals]
    
    logist = calc_logist(y_offset)

    (a, b), covariance = curve_fit(logist,
                                   x,
                                   quant_vals,
                                   bounds=((2*yint_min, 0),
                                           (2*yint_max, np.inf)))
    return a, b


def get_slope(fit):
    return fit[0]


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args: {}".format(args))

    contin_renorm_main(args)
