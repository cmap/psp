import numpy as np
import gct

def update_prov_code(new_entry, existing_prov_code):
    pass


def filter_samples(df, sample_pct_cutoff):
    # Return df (potentially of different size than original df)
    pass


def filter_probes(df, probe_pct_cutoff, probe_sd_cutoff):
    # Identify rows manually labeled for rejection
    # Return df (potentially of different size than original df)
    pass

def optimize_sample_balance(blah):
    pass

def main():
    gct_obj = gct.GCT("/Users/lev/code/PSP/python/functional_tests/test_p100.gct")
    gct_obj.read(row_inds=range(100),col_inds=range(10))
    data = gct_obj.matrix
    print data


# if __name__ == '__main__':
#     args = build_parser().parse_args(sys.argv[1:])
#     setup_logger.setup(args.verbose, args.log_file)
#     logger.debug("args:  {}".format(args))
#
#     make_args_abs_paths(args)
#     logger.debug("args after make_args_abs_paths args:  {}".format(args))
#
#     main(args)