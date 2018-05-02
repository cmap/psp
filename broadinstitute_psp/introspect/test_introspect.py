import logging
import pandas as pd
import os
import unittest
import cmapPy.pandasGEXpress.parse as parse
import cmapPy.pandasGEXpress.setup_GCToo_logger as setup_logger
import broadinstitute_psp.introspect.introspect as introspect

FUNCTIONAL_TESTS_DIR = "introspect/functional_tests"

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestIntrospect(unittest.TestCase):

	def test_main1(self):
		input_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR,
		                              "test_introspect_main.gct")
		output_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR,
		                               "test_introspect_main_out.gct")
		expected_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR,
		                                 "test_introspect_main_expected.gct")

		args_string = "-i {} -o {} -fa chd1".format(input_gct_path, output_gct_path)
		args = introspect.build_parser().parse_args(args_string.split())

		introspect.main(args)

		# Read in output and expected gcts and confirm that they're equal
		output_gct = parse.parse(output_gct_path)
		expected_gct = parse.parse(expected_gct_path)

		pd.util.testing.assert_almost_equal(expected_gct.data_df, output_gct.data_df, check_less_precise=2)
		pd.testing.assert_frame_equal(expected_gct.row_metadata_df, output_gct.row_metadata_df)
		pd.testing.assert_frame_equal(expected_gct.col_metadata_df, output_gct.col_metadata_df)

		# Clean up
		os.remove(output_gct_path)

	def test_main2(self):
		input_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR,
		                              "test_introspect_main.gct")
		output_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR,
		                               "test_introspect_main_out2.gct")
		expected_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR,
		                                 "test_introspect_main_expected2.gct")

		args_string = "-i {} -o {} -fa moa".format(
			input_gct_path, output_gct_path)
		args = introspect.build_parser().parse_args(args_string.split())

		introspect.main(args)

		# Read in output and expected gcts and confirm that they're equal
		output_gct = parse.parse(output_gct_path)
		expected_gct = parse.parse(expected_gct_path)

		pd.util.testing.assert_almost_equal(expected_gct.data_df, output_gct.data_df, check_less_precise=True)
		pd.util.testing.assert_frame_equal(expected_gct.row_metadata_df, output_gct.row_metadata_df)
		pd.util.testing.assert_frame_equal(expected_gct.col_metadata_df, output_gct.col_metadata_df)

		# Clean up
		os.remove(output_gct_path)


if __name__ == '__main__':
	setup_logger.setup(verbose=True)
	unittest.main()

