"""
introspect.py

Compute all within-dataset connectivities.

"""

import logging
import sys
import argparse

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.write_gct as wg
import broadinstitute_psp.steep.steep as steep
import broadinstitute_psp.sip.sip as sip

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
	parser.add_argument("--in_gct_path", "-i", required=True,
						help="path to gct file of profiles")

	# Optional args
	parser.add_argument("--out_sip_name", "-osi", default="sip_output.gct",
						help="what to name the output connectivity file")
	parser.add_argument("--similarity_metric", "-s", default="spearman",
						choices=["spearman", "pearson"],
						help="similarity metric to use for comparing columns")
	parser.add_argument("--connectivity_metric", "-c", default="ks_test",
						choices=["ks_test", "percentile_score"],
						help="metric to use for computing connectivity")
	parser.add_argument("--fields_to_aggregate", "-fa",
						nargs="+", default=["pert_id", "cell_id", "pert_time"],
						help="list of metadata fields to use in aggregating replicates")
	parser.add_argument("--verbose", "-v", action="store_true", default=False,
						help="whether to increase the # of messages reported")

	return parser


def main(args):

	gct = parse.parse(args.in_gct_path)

	(_, conn_gct) = do_steep_and_sip(
		gct, args.similarity_metric,
		args.connectivity_metric, args.fields_to_aggregate)

	# Write output gct
	wg.write(conn_gct, args.out_sip_name, data_null="NaN", filler_null="NaN", metadata_null="NaN")


def do_steep_and_sip(gct, similarity_metric, connectivity_metric, fields_to_aggregate):
	""" Perform steep and sip on the same GCT. AKA introspect.

	Args:
	    gct:
	    similarity_metric:
	    connectivity_metric:
	    fields_to_aggregate:

	Returns:
	    sim_gct
	    conn_gct

	"""

	#----------STEEP--------#

	sim_df = steep.compute_similarity_within_df(gct.data_df, similarity_metric)

	# Row and column metadata are both from gct
	metadata_df = gct.col_metadata_df

	# Append column to metadata_df indicating which similarity_metric was used
	metadata_df[SIMILARITY_METRIC_FIELD] = similarity_metric

	# Assemble similarity gct
	sim_gct = GCToo.GCToo(data_df=sim_df, row_metadata_df=metadata_df,
	                      col_metadata_df=metadata_df)

	#----------SIP----------#

	# Check symmetry
	(is_test_df_sym, _) = sip.check_symmetry(sim_gct.data_df, sim_gct.data_df)

	# Create deep copies of sim_gct in order to leave the original GCT untouched
	test_gct = GCToo.GCToo(data_df=sim_df.copy(deep=True),
	                       row_metadata_df=metadata_df.copy(deep=True),
	                       col_metadata_df=metadata_df.copy(deep=True))
	bg_gct = GCToo.GCToo(data_df=sim_df.copy(deep=True),
	                     row_metadata_df=metadata_df.copy(deep=True),
	                     col_metadata_df=metadata_df.copy(deep=True))

	# Create an aggregated metadata field for index and columns of sim_gct
	# and sort by that field
	(test_gct, bg_gct) = sip.create_aggregated_fields_in_GCTs(
		test_gct, bg_gct, fields_to_aggregate, fields_to_aggregate,
		fields_to_aggregate, QUERY_FIELD_NAME, TARGET_FIELD_NAME, SEPARATOR)

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