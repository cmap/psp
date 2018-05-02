import logging
import unittest
import os
import pandas as pd

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.parse as parse
import annotate_gct_from_mapping as agfm

functional_tests_dir = "utils/functional_tests/"

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestAnnotate(unittest.TestCase):

    def test_main(self):
        gct_path = os.path.join(functional_tests_dir, "test_annotate_gct_from_mapping_in.gct")
        mapping_path = os.path.join(functional_tests_dir, "test_annotate_gct_from_mapping.tsv")
        expected_gct_path = os.path.join(functional_tests_dir, "test_annotate_gct_from_mapping_expected.gct")
        out_path = os.path.join(functional_tests_dir, "test_annotate_gct_from_mapping_out.gct")

        args_string = "-i {} -m {} -o {} -f {}".format(
            gct_path, mapping_path, out_path, "pert_iname")
        args = agfm.build_parser().parse_args(args_string.split())

        agfm.main(args)

        # Read in expected and actual outputs
        e_gct = parse.parse(expected_gct_path)
        out_gct = parse.parse(out_path)

        pd.util.testing.assert_frame_equal(e_gct.data_df, out_gct.data_df)
        pd.util.testing.assert_frame_equal(e_gct.row_metadata_df, out_gct.row_metadata_df)
        pd.util.testing.assert_frame_equal(e_gct.col_metadata_df, out_gct.col_metadata_df)

        # Clean up
        os.remove(out_path)

    def test_annotate_meta_df(self):
        meta_df = pd.DataFrame(
            [["a"], ["b"], ["c"]], index=["A", "B", "C"], columns=["pert_iname"])

        # Make a copy of meta_df to use for second test
        meta_df2 = meta_df.copy(deep=True)

        to_entries = pd.Series(
            ["inhibitor", "activator", "killer", "winner"],
            index=["f", "b", "c", "d"])
        to_entries.name = "moa"

        e_meta_df = pd.DataFrame(
            [["a", "NA"], ["b", "activator"], ["c", "killer"]],
            index=["A", "B", "C"], columns=["pert_iname", "moa"])

        agfm.annotate_meta_df(meta_df, to_entries, "pert_iname", "NA")
        pd.util.testing.assert_frame_equal(meta_df, e_meta_df)

        with self.assertRaises(AssertionError) as e:
            agfm.annotate_meta_df(meta_df, to_entries, "thing3", "NA")
        self.assertIn("gct_from_field must be a metadata header", str(e.exception))

        different_to_entries = pd.Series(
            ["X", "Y", "Z"],
            index=["A", "C", "B"])
        different_to_entries.name = "letter"

        e_meta_df2 = pd.DataFrame(
            [["a", "X"], ["b", "Z"], ["c", "Y"]],
            index=["A", "B", "C"], columns=["pert_iname", "letter"])

        agfm.annotate_meta_df(meta_df2, different_to_entries, None, "NA")
        pd.util.testing.assert_frame_equal(meta_df2, e_meta_df2)


if __name__ == '__main__':
    setup_logger.setup(verbose=True)
    unittest.main()
