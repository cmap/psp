"""
steep_and_sip_external_on_clue.py

Most of the arguments come from the config file.

Runs steep_and_sip_external for each cell line in the corpus.

Requirements:
Run steep_and_sip_external.py on all cell lines.
Parse external_gct_path only once.
Create appropriate error message.
Write text file indicating success or failure.


"""

import argparse
import ConfigParser
import datetime
import logging
import os
import sys
import time

import broadinstitute_cmap.io.pandasGEXpress.parse as pg
import broadinstitute_cmap.io.pandasGEXpress.write_gct as wg
import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.sip.steep_and_sip_external as steep_and_sip_external

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)

INTERNAL_GCT_FORMAT = "{assay}_{cell}_QCNORM.gct"
BG_GCT_FORMAT = "{assay}_{cell}_SIM.gct"
OUT_STEEP_FORMAT = "{cell}_SIM.gct"
OUT_SIP_FORMAT = "{cell}_CONN.gct"


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--assay", "-a", required=True, choices=["GCP", "P100"],
                        help="which assay's data to query")
    parser.add_argument("--external_gct_path", "-e", required=True,
                        help="path to gct file of external profiles")
    parser.add_argument("--out_dir", "-o", required=True,
                        help="output directory")
    parser.add_argument("--psp_on_clue_config_path", "-p",
                        required=True,
                        help="filepath to psp_on_clue.cfg")

    # Optional args
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):

    try:

        # Read and unpack config file
        (cells, internal_gct_dir, bg_gct_dir, fields_to_aggregate_for_external_profiles,
         fields_to_aggregate_for_internal_profiles, similarity_metric,
         connectivity_metric) = read_config_file(args.psp_on_clue_config_path)

        # Read in the external profiles only once
        external_gct = pg.parse(args.external_gct_path, convert_neg_666=False, make_multiindex=True)

        # Loop over cell lines in corpus
        for cell in cells:

            # Import gct with the internal profiles for this cell line
            internal_gct_path = os.path.join(internal_gct_dir, INTERNAL_GCT_FORMAT.format(
                args.assay, cell))
            internal_gct = pg.parse(internal_gct_path, convert_neg_666=False, make_multiindex=True)

            # Import gct with the similarity matrix for this cell line
            bg_gct_path = os.path.join(bg_gct_dir, BG_GCT_FORMAT.format(
                args.assay, cell))
            bg_gct = pg.parse(bg_gct_path, convert_neg_666=False, make_multiindex=True)

            try:
                (sim_gct, conn_gct) = steep_and_sip_external.do_steep_and_sip(
                    external_gct, internal_gct, bg_gct, "spearman",
                    "ks_test", fields_to_aggregate_for_external_profiles,
                    fields_to_aggregate_for_internal_profiles)

                # Write output gcts
                out_steep_name = os.path.join(args.out_dir, OUT_STEEP_FORMAT.format(cell))
                out_sip_name = os.path.join(args.out_dir, OUT_SIP_FORMAT.format(cell))
                wg.write(sim_gct, out_steep_name, data_null="NaN", metadata_null="NA", filler_null="NA")
                wg.write(conn_gct, out_sip_name, data_null="NaN", filler_null="NaN", metadata_null="NaN")

            except Exception as e:
                msg = "Could not complete for cell {}. stacktrace: {}".format(cell, str(e))

                # Write failure.txt
                write_text_file(os.path.join(args.out_dir, "failure.txt"), msg)

                # Raise exception
                logger.error(msg)
                raise Exception(msg)

        # Write success.txt with timestamp
        timestamp = datetime.datetime.fromtimestamp(time.time()).strftime('%Y-%m-%d %H:%M:%S')
        write_text_file(os.path.join(args.out_dir, "success.txt"), timestamp)

    except Exception as e:
        msg = "Something failed outside of looping over cell lines. stacktrace: "

        # Write failure.txt
        write_text_file(os.path.join(args.out_dir, "failure.txt"), msg)

        # Raise exception
        logger.error(msg)
        raise Exception(msg)


def read_config_file(config_path):

    assert os.path.exists(config_path), (
        "Config file can't be found. config_path: {}".format(config_path))

    # Read config file
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_path)

    # Return config fields as dictionarires
    config_corpus = dict(config_parser.items("corpus"))
    config_metadata = dict(config_parser.items("metadata"))
    config_algorithms = dict(config_parser.items("algorithms"))

    # Unpack the config file
    cells = eval(config_corpus["cells"])
    internal_gct_dir = config_corpus["qcnorm_dir"]
    bg_gct_dir = config_corpus["sim_dir"]
    fields_to_aggregate_for_external_profiles = eval(config_metadata["fields_to_aggregate_for_external_profiles"])
    fields_to_aggregate_for_internal_profiles = eval(config_metadata["fields_to_aggregate_for_internal_profiles"])
    similarity_metric = config_algorithms["similarity_metric"]
    connectivity_metric = config_algorithms["connectivity_metric"]

    return cells, internal_gct_dir, bg_gct_dir, \
           fields_to_aggregate_for_external_profiles, \
           fields_to_aggregate_for_internal_profiles, \
           similarity_metric, connectivity_metric


def write_text_file(file_name, file_contents):
    f = open(file_name, 'w')
    f.write(file_contents)
    f.close()


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)