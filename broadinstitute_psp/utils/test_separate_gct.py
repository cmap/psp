import logging
import os
import pandas as pd
import unittest

import cmapPy.pandasGEXpress.parse as parse
import broadinstitute_psp.utils.separate_gct as sg
import broadinstitute_psp.utils.setup_logger as setup_logger

logger = logging.getLogger(setup_logger.LOGGER_NAME)

functional_tests_dir = "utils/functional_tests/"

in_gct_path = functional_tests_dir + "test_separate_in.gct"
thing1_gct_path = functional_tests_dir + "test_separate_expected_thing1.gct"
thing2_gct_path = functional_tests_dir + "test_separate_expected_thing2.gct"
a375_gct_path = functional_tests_dir + "test_separate_expected_A375.gct"
ht29_gct_path = functional_tests_dir + "test_separate_expected_HT29.gct"
a549_gct_path = functional_tests_dir + "test_separate_expected_A549.gct"

in_gct = parse.parse(in_gct_path)
thing1_gct = parse.parse(thing1_gct_path)
thing2_gct = parse.parse(thing2_gct_path)
a375_gct = parse.parse(a375_gct_path)
ht29_gct = parse.parse(ht29_gct_path)
a549_gct = parse.parse(a549_gct_path)


class TestSeparateGct(unittest.TestCase):
    def test_separate_row(self):
        (thing_gcts, thing_fields) = sg.separate(in_gct, "thing", "row")

        self.assertListEqual(thing_fields, [1, 2])

        pd.util.testing.assert_frame_equal(thing_gcts[0].data_df, thing1_gct.data_df)
        pd.util.testing.assert_frame_equal(thing_gcts[0].row_metadata_df, thing1_gct.row_metadata_df)
        pd.util.testing.assert_frame_equal(thing_gcts[0].col_metadata_df, thing1_gct.col_metadata_df)

        pd.util.testing.assert_frame_equal(thing_gcts[1].data_df, thing2_gct.data_df)
        pd.util.testing.assert_frame_equal(thing_gcts[1].row_metadata_df, thing2_gct.row_metadata_df)
        pd.util.testing.assert_frame_equal(thing_gcts[1].col_metadata_df, thing2_gct.col_metadata_df)

    def test_separate_col(self):
        (cell_id_gcts, cell_id_fields) = sg.separate(in_gct, "cell_id", "col")

        self.assertListEqual(cell_id_fields, ["A375", "HT29", "A549"])

        pd.util.testing.assert_frame_equal(cell_id_gcts[0].data_df, a375_gct.data_df)
        pd.util.testing.assert_frame_equal(cell_id_gcts[0].row_metadata_df, a375_gct.row_metadata_df)
        pd.util.testing.assert_frame_equal(cell_id_gcts[0].col_metadata_df, a375_gct.col_metadata_df)

        pd.util.testing.assert_frame_equal(cell_id_gcts[1].data_df, ht29_gct.data_df)
        pd.util.testing.assert_frame_equal(cell_id_gcts[1].row_metadata_df, ht29_gct.row_metadata_df)
        pd.util.testing.assert_frame_equal(cell_id_gcts[1].col_metadata_df, ht29_gct.col_metadata_df)

        pd.util.testing.assert_frame_equal(cell_id_gcts[2].data_df, a549_gct.data_df)
        pd.util.testing.assert_frame_equal(cell_id_gcts[2].row_metadata_df, a549_gct.row_metadata_df)
        pd.util.testing.assert_frame_equal(cell_id_gcts[2].col_metadata_df, a549_gct.col_metadata_df)

    def test_main_row(self):
        out_prefix = "test_main_row_"
        args_string = "-i {} -sf {} -rc {} -od {} -op {}".format(
            in_gct_path, "thing", "row", functional_tests_dir, out_prefix)
        args = sg.build_parser().parse_args(args_string.split())

        # Run main method
        sg.main(args)

        # Clean up
        os.remove(os.path.join(functional_tests_dir, out_prefix + "1.gct"))
        os.remove(os.path.join(functional_tests_dir, out_prefix + "2.gct"))

    def test_main_col(self):
        out_prefix = "test_main_col_"
        args_string = "-i {} -sf {} -rc {} -od {} -op {}".format(
            in_gct_path, "cell_id", "col", functional_tests_dir, out_prefix)
        args = sg.build_parser().parse_args(args_string.split())

        # Run main method
        sg.main(args)

        # Clean up
        os.remove(os.path.join(functional_tests_dir, out_prefix + "A375.gct"))
        os.remove(os.path.join(functional_tests_dir, out_prefix + "HT29.gct"))
        os.remove(os.path.join(functional_tests_dir, out_prefix + "A549.gct"))


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
