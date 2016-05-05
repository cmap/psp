import logging
import utils.setup_logger as setup_logger
import sys
import argparse
import in_out.parse_gctoo as parse_gctoo
import pandas as pd

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

"""
Converts a QC gct file into a .pw (plate-well) file that can be used
easily with Plato. Returns median, mean, MAD, and sd of the values in
each column (or sample) of the QC file.

Input is a gct file. Output is a pw file.

"""

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("gct_file_path", type=str, help="filepath to gct file")
    parser.add_argument("out_pw_file_path", type=str,
                        help="filepath to output pw file")
    
    # Optional args
    parser.add_argument("-plate_field", type=str, default="det_plate",
                        help="metadata field name specifying the plate")
    parser.add_argument("-well_field", type=str, default="det_well",
                        help="metadata field name specifying the well")

    return parser


def main(args):
    # Import gct
    gct = parse_gctoo.parse(args.gct_file_path)

    # Get plate and well names
    (plate_names, well_names) = extract_plate_and_well_names(
        gct.col_metadata_df, args.plate_field, args.well_field)

    # TODO(lev): make this its own method

    # Calculate metrics for each sample
    medium_over_heavy_medians = gct.data_df.median(axis=0).values
    medium_over_heavy_means = gct.data_df.mean(axis=0).values
    medium_over_heavy_mads = gct.data_df.mad(axis=0).values
    medium_over_heavy_sds = gct.data_df.std(axis=0).values

    # Assemble plate_names, well_names, and metrics into a dataframe
    temp_df = pd.DataFrame.from_dict({
        "plate_name": plate_names,
        "well_name": well_names,
        "medium_over_heavy_median": medium_over_heavy_medians,
        "medium_over_heavy_mean": medium_over_heavy_means,
        "medium_over_heavy_mad": medium_over_heavy_mads,
        "medium_over_heavy_sd": medium_over_heavy_sds
    })

    # Specify order of columns
    cols = ["plate_name", "well_name", "medium_over_heavy_median",
            "medium_over_heavy_mean", "medium_over_heavy_mad",
            "medium_over_heavy_sd"]
    out_df = temp_df[cols]

    # Write to pw file
    out_df.to_csv(args.out_pw_file_path, sep="\t", na_rep="NaN", index=False)
    logger.info("PW file written to {}".format(args.out_pw_file_path))


def extract_plate_and_well_names(col_meta, plate_field, well_field):

    # Extract plate metadata
    plate_names = col_meta[plate_field].values

    # Verify that all samples have the same plate name
    plate_names_same = True
    for plate in plate_names:
        plate_names_same = (plate_names_same and plate == plate_names[0])
        assert plate_names_same

    # Extract well metadata
    well_names = col_meta[well_field].values

    return plate_names, well_names


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup()
    main(args)
