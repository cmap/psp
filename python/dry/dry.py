import numpy as np
import pandas as pd
import in_out.gct as gct


def determine_assay_type(prov_code):
    """
    Convert assay_type from string to numpy array and compare to a list of allowable entries.

    Input
        -prov_code: string
    Output
        -assay_type: numpy string array, length=1
    """
    allowed_assay_types = np.array(["PR1", "GR1", "DR1"], dtype=np.str)
    assay_type = np.array(prov_code, dtype=np.str)

    # Check that assay type is one we expect
    if assay_type not in allowed_assay_types:
        error_msg = ("The first entry of the provenance code should be the assay type. " +
                     "Assay type not recognized. assay_type: {}")
        raise Exception(error_msg.format(assay_type))

    return assay_type


def update_prov_code(new_entry, existing_prov_code):
    """
    Append entry to provenance code.

    Input
        -new_entry: string
        -existing_prov_code: numpy string array
    Output
        -updated_prov_code : numpy string array
    """
    updated_prov_code = np.append(existing_prov_code, new_entry)
    return updated_prov_code

def parse_prov_code(prov_code):
    # Not sure when I'll use this, but seems good to keep
    # prov_code = "PR1+L2X+SF8+PF8.5"
    # expected_out = np.array(["PR1", "L2X", "SF8", "PF8.5"], dtype=np.str)
    # actual_out = dry.parse_prov_code(prov_code)
    # self.assertTrue(np.array_equal(actual_out, expected_out),
    #                 "Expected output: {}, Actual output: {}".format(expected_out, actual_out))
    # # Split along + in case there are other entries in the prov_code already
    # prov_code_chunks = np.array(prov_code.split("+"))
    pass


def filter_samples(df, sample_pct_cutoff):
    """
    Remove samples (i.e. columns) with more than sample_pct_cutoff NaN values.

    Input
        -df: pandas dataframe
        -sample_pct_cutoff: number
    Output
        -out_df: pandas dataframe (potentially smaller than original df)
    """
    # Number of NaNs per sample
    num_nans = df.isnull().sum()

    # Number of probes
    num_probes = df.shape[1]

    # Percent NaNs per sample
    pct_nans_per_sample = num_nans.divide(num_probes)

    # Only return samples with fewer % of NaNs than the cutoff
    out_df = df.loc[:, pct_nans_per_sample < sample_pct_cutoff]
    return out_df.values


def filter_probes(df, probe_pct_cutoff, probe_sd_cutoff):
    """
    Remove probes (i.e. rows) with more than probe_pct_cutoff NaN values.
    Also remove probes with standard deviation higher than probe_sd_cutoff.

    Input
        -df: pandas dataframe
        -probe_pct_cutoff: number
    Output
        -out_df: pandas dataframe (potentially smaller than original df)
    """

    ### CONTINUE HERE. ###

    # Number of NaNs per sample
    num_nans = df.isnull().sum()

    # Number of probes
    num_probes = df.shape[1]

    # Percent NaNs per sample
    pct_nans_per_sample = num_nans.divide(num_probes)

    # Only return samples with fewer % of NaNs than the cutoff
    out_df = df.loc[:, pct_nans_per_sample < probe_pct_cutoff]
    return out_df


    # num_nans_other_way = df.isnull().sum(axis=1)
    # print num_nans_other_way
    #     # Number of NaNs per sample
    # num_nans = df.isnull().sum()
    #
    # # Number of probes
    # num_probes = df.shape[1]
    #
    # # Percent NaNs per sample
    # pct_nans_per_sample = num_nans.divide(num_probes)
    #
    # # Only return samples with fewer % of NaNs than the cutoff
    # out_df = df.loc[:, pct_nans_per_sample < sample_pct_cutoff]
    # return out_df
    #
    #
    # # Identify rows manually labeled for rejection
    # # Return df (potentially of different size than original df)
    # pass

def manual_probe_rejection(df):
    """
    Also remove those probes that were manually selected for rejection.
    :param df:
    :return:
    """



def optimize_sample_balance(blah):
    pass


def main():
    # gct_obj = gct.GCT("/Users/lev/code/PSP/python/functional_tests/test_p100.gct")
    # print(gct_obj.version)
    # gct_obj.read(row_inds=range(100),col_inds=range(10))
    # data = gct_obj.matrix
    # print(data)
    pass


# if __name__ == '__main__':
#     args = build_parser().parse_args(sys.argv[1:])
#     setup_logger.setup(args.verbose, args.log_file)
#     logger.debug("args:  {}".format(args))
#
#     make_args_abs_paths(args)
#     logger.debug("args after make_args_abs_paths args:  {}".format(args))
#
#     main(args)