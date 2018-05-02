"""
external_query_many.py

Runs external_query.py for each cell line in the corpus and optionally also
computes introspect. Most of the arguments come from the config file.
The default config file points to the latest signature and similarity
directories.

"""

import ConfigParser
import argparse
import datetime
import logging
import os
import sys
import traceback

import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.write_gct as wg
import cmapPy.pandasGEXpress.concat as cg

import broadinstitute_psp.external_query.external_query as eq
import broadinstitute_psp.introspect.introspect as introspect
import broadinstitute_psp.utils.setup_logger as setup_logger

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)

INTERNAL_GCT_FORMAT = "{assay}_{cell}_DIFF.gct"
BG_GCT_FORMAT = "{assay}_{cell}_SIM.gct"
OUT_STEEP_FORMAT = "{cell}_SIM.gct"
OUT_SIP_FORMAT = "{cell}_CONN.gct"
OUT_CONCATED_NAME = "CONCATED_CONN.gct"
OUT_INTROSPECT_NAME = "INTROSPECT_CONN.gct"


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--assay", "-a", required=True, choices=["GCP", "P100"],
                        help="which assay's data to query")
    parser.add_argument("--introspect", "-i", action="store_true", default=False,
                        help="whether to also compute introspect")
    parser.add_argument("--external_gct_path", "-e", required=True,
                        help="path to gct file of external profiles")
    parser.add_argument("--out_dir", "-o", required=True,
                        help="directory in which to dump output")
    parser.add_argument("--psp_on_clue_config_path", "-p",
                        default="clue/psp_on_clue.yml",
                        help="filepath to psp_on_clue.yml")
    parser.add_argument("--fields_to_aggregate_for_external_profiles", "-fae",
                        nargs="*", default=["pert_id", "cell_id", "pert_time"],
                        help="list of metadata fields to use in aggregating replicates in external profiles")
    parser.add_argument("--all", default=False, action="store_true",
                        help="whether to produce all output matrices")

    # Optional args
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):

    # Record start_time
    start_time = datetime.datetime.now()
    start_time_msg = "external_query_many.py started at {}".format(
        start_time.strftime('%Y-%m-%d %H:%M:%S'))

    # Create output directory
    assert os.path.exists(args.out_dir), "args.out_dir: {}".format(args.out_dir)

    try:

        # Read and unpack config file
        (cells, internal_gct_dir, bg_gct_dir,
         fields_to_aggregate_for_internal_profiles, similarity_metric,
         connectivity_metric) = read_config_file(args.psp_on_clue_config_path)

        # Read in the external profiles only once
        external_gct = parse.parse(args.external_gct_path)

        # If requested, do introspect
        (_, introspect_gct) = introspect.do_steep_and_sip(
            external_gct, similarity_metric, connectivity_metric,
            args.fields_to_aggregate_for_external_profiles)

        # Write introspect result
        actual_out_introspect_name = os.path.join(args.out_dir, OUT_INTROSPECT_NAME)
        wg.write(introspect_gct, actual_out_introspect_name, data_null="NaN", metadata_null="NaN", filler_null="NaN")

        # Initialize list to store connectivity gcts
        list_of_conn_gcts = []

        # Loop over cell lines in corpus
        for cell in cells:

            # Import gct with the internal profiles for this cell line
            internal_gct_path = os.path.join(internal_gct_dir, INTERNAL_GCT_FORMAT.format(
                assay=args.assay, cell=cell))
            internal_gct = parse.parse(internal_gct_path)

            # Import gct with the similarity matrix for this cell line
            bg_gct_path = os.path.join(bg_gct_dir, BG_GCT_FORMAT.format(
                assay=args.assay, cell=cell))
            bg_gct = parse.parse(bg_gct_path)

            (sim_gct, conn_gct) = eq.do_steep_and_sip(
                external_gct, internal_gct, bg_gct, "spearman",
                "ks_test", args.fields_to_aggregate_for_external_profiles,
                fields_to_aggregate_for_internal_profiles)

            # Append this connectivity gct
            list_of_conn_gcts.append(conn_gct)

            # Write all output gcts if requested
            if args.all:
                out_steep_name = os.path.join(args.out_dir, OUT_STEEP_FORMAT.format(cell=cell))
                out_sip_name = os.path.join(args.out_dir, OUT_SIP_FORMAT.format(cell=cell))

                wg.write(sim_gct, out_steep_name)
                wg.write(conn_gct, out_sip_name)

        # Concatenate connectivity GCTs
        concated = cg.vstack(list_of_conn_gcts)
        actual_out_concated_name = os.path.join(args.out_dir, OUT_CONCATED_NAME)

        # Write concatenated result
        wg.write(concated, actual_out_concated_name, data_null="NaN", filler_null="NaN", metadata_null="NaN")

        # Write success.txt with timestamp
        success_path = os.path.join(args.out_dir, "success.txt")
        write_success(success_path, start_time_msg)

        # Return how much time it took
        end_time = datetime.datetime.now()
        seconds_elapsed = (end_time - start_time).seconds
        logger.info("external_query_many.py completed in {:.0f} sec.".format(seconds_elapsed))

    except Exception:
        failure_path = os.path.join(args.out_dir, "failure.txt")
        msg = "external_query_many.py failed. See {} for stacktrace.".format(failure_path)

        # Write failure.txt
        write_failure(failure_path, start_time_msg)

        # Raise exception
        logger.error(msg)
        raise Exception(msg)

    return None


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
    internal_gct_dir = config_corpus["signature_dir"]
    bg_gct_dir = config_corpus["sim_dir"]
    fields_to_aggregate_for_internal_profiles = eval(config_metadata["fields_to_aggregate_for_internal_profiles"])
    similarity_metric = config_algorithms["similarity_metric"]
    connectivity_metric = config_algorithms["connectivity_metric"]

    return cells, internal_gct_dir, bg_gct_dir, \
           fields_to_aggregate_for_internal_profiles, \
           similarity_metric, connectivity_metric


def write_success(file_name, start_time_msg):
    # Create timestamp
    timestamp = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')

    # Write timestamp to file_name
    f = open(file_name, 'w')
    f.write(start_time_msg + "\n")
    f.write("external_query_many.py completed at {}\n".format(timestamp))
    f.close()


def write_failure(file_name, start_time_msg):
    # Record stacktrace
    _, _, exc_traceback = sys.exc_info()

    # Write stacktrace to file_name
    f = open(file_name, "w")
    f.write(start_time_msg + "\n")
    traceback.print_exc(exc_traceback, file=f)
    f.close()


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)
