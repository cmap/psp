"""
external_query.py

Run steep.py (compute similarities) and sip.py (compute connectivities) on an
external gct file. The required inputs are a path to a gct of external profiles
(probes x samples), a path to a gct of internal profiles, and a path to a gct
of the pre-computed similarity matrix of the internal profiles against
themselves.

N.B. The internal_gcts and bg_gcts should each contain only WITHIN-cell
calculations. In other words, this script needs to be run 6 times to get
results for each cell line in the corpus.

"""

import logging
import sys
import argparse

import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.steep.steep as steep
import broadinstitute_psp.sip.sip as sip
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.write_gct as wg

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)

SIMILARITY_METRIC_FIELD = "similarity_metric"
CONNECTIVITY_METRIC_FIELD = "connectivity_metric"
QUERY_FIELD_NAME = "query_field"
TARGET_FIELD_NAME = "target_field"
SEPARATOR = ":"


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--external_gct_path", "-e", required=True,
                        help="path to gct file of external profiles")
    parser.add_argument("--internal_gct_path", "-i", required=True,
                        help="path to gct file of internal profiles")
    parser.add_argument("--bg_gct_path", "-b", required=True,
                        help="path to background similarity gct file")

    # Optional args
    parser.add_argument("--out_steep_name", "-ost", default="steep_output.gct",
                        help="what to name the output similarity file")
    parser.add_argument("--out_sip_name", "-osi", default="sip_output.gct",
                        help="what to name the output connectivity file")
    parser.add_argument("--similarity_metric", "-s", default="spearman",
                        choices=["spearman", "pearson"],
                        help="metric to use for comparing sample profiles")
    parser.add_argument("--connectivity_metric", "-c", default="ks_test",
                        choices=["ks_test", "percentile_score"],
                        help="metric to use for computing connectivity")
    parser.add_argument("--fields_to_aggregate_for_external_profiles", "-fae",
                        nargs="+", default=["pert_id", "cell_id", "pert_time"],
                        help="list of metadata fields to use in aggregating replicates in external profiles")
    parser.add_argument("--fields_to_aggregate_for_internal_profiles", "-fai",
                        nargs="+", default=["pert_id", "cell_id", "pert_time"],
                        help="list of metadata fields to use in aggregating replicates in internal profiles")
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    return parser


def main(args):

    # Parse input gcts
    external_gct = parse.parse(args.external_gct_path)
    internal_gct = parse.parse(args.internal_gct_path)
    bg_gct = parse.parse(args.bg_gct_path)

    # Meat of the script
    (sim_gct, conn_gct) = do_steep_and_sip(
        external_gct, internal_gct, bg_gct, args.similarity_metric,
        args.connectivity_metric,
        args.fields_to_aggregate_for_external_profiles,
        args.fields_to_aggregate_for_internal_profiles)

    # Write output gcts
    wg.write(sim_gct, args.out_steep_name, data_null="NaN", metadata_null="NaN", filler_null="NaN")
    wg.write(conn_gct, args.out_sip_name, data_null="NaN", filler_null="NaN", metadata_null="NaN")


def do_steep_and_sip(external_gct, internal_gct, bg_gct, similarity_metric,
                     connectivity_metric,
                     fields_to_aggregate_for_external_profiles,
                     fields_to_aggregate_for_internal_profiles):

    #----------STEEP----------#

    # Compute similarity between external and internal profiles
    sim_df = steep.compute_similarity_bw_two_dfs(internal_gct.data_df,
                                                 external_gct.data_df,
                                                 similarity_metric)

    # Row metadata is from gct1, column metadata is from gct2
    row_metadata_for_sim_df = internal_gct.col_metadata_df
    col_metadata_for_sim_df = external_gct.col_metadata_df

    # Append column to both metadata_dfs indicating which similarity_metric was used
    row_metadata_for_sim_df[SIMILARITY_METRIC_FIELD] = similarity_metric
    col_metadata_for_sim_df[SIMILARITY_METRIC_FIELD] = similarity_metric

    # Assemble similarity gct
    sim_gct = GCToo.GCToo(sim_df, row_metadata_for_sim_df, col_metadata_for_sim_df)

    #----------SIP----------#

    # Check symmetry
    (is_test_df_sym, is_bg_df_sym) = sip.check_symmetry(sim_gct.data_df, bg_gct.data_df)

    # Create an aggregated metadata field for index and columns of both gcts
    # and sort by that field
    (test_gct, bg_gct) = sip.create_aggregated_fields_in_GCTs(
        sim_gct, bg_gct,
        fields_to_aggregate_for_external_profiles,
        fields_to_aggregate_for_internal_profiles,
        fields_to_aggregate_for_internal_profiles,
        QUERY_FIELD_NAME, TARGET_FIELD_NAME, SEPARATOR)

    # Compute connectivity
    (_, signed_conn_gct) = sip.compute_connectivities(
        test_gct, bg_gct, QUERY_FIELD_NAME, TARGET_FIELD_NAME, TARGET_FIELD_NAME,
        connectivity_metric, is_test_df_sym, SEPARATOR)

    # Append to queries a new column saying what connectivity metric was used
    sip.add_connectivity_metric_to_metadata(signed_conn_gct.col_metadata_df, connectivity_metric, CONNECTIVITY_METRIC_FIELD)
    sip.add_connectivity_metric_to_metadata(signed_conn_gct.row_metadata_df, connectivity_metric, CONNECTIVITY_METRIC_FIELD)

    return sim_gct, signed_conn_gct


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)