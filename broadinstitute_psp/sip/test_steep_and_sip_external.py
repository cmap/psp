import unittest
import os
import pandas as pd

import broadinstitute_cmap.io.pandasGEXpress.parse_gct as pg
import steep_and_sip_external as sse

FUNCTIONAL_TESTS_DIR = "sip/functional_tests"


class TestSteepAndSipExternal(unittest.TestCase):
    def test_main(self):
        external_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_external.gct")
        internal_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_internal.gct")
        bg_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_bg.gct")
        expected_sim_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_expected_sim.gct")
        expected_conn_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_expected_conn.gct")

        out_steep_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_steep_output.gct")
        out_sip_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_sip_output.gct")

        args_string = "-e {} -i {} -b {} -ost {} -osi {} -p {}".format(
            external_gct_path, internal_gct_path, bg_gct_path,
            out_steep_path, out_sip_path, "psp_production.cfg"
        )
        args = sse.build_parser().parse_args(args_string.split())
        sse.main(args)

        expected_sim_gct = pg.parse(expected_sim_path)
        expected_conn_gct = pg.parse(expected_conn_path)
        actual_sim_gct = pg.parse(out_steep_path)
        actual_conn_gct = pg.parse(out_sip_path)

        pd.util.testing.assert_frame_equal(expected_sim_gct.data_df, actual_sim_gct.data_df)
        pd.util.testing.assert_frame_equal(expected_sim_gct.row_metadata_df, actual_sim_gct.row_metadata_df)
        pd.util.testing.assert_frame_equal(expected_sim_gct.col_metadata_df, actual_sim_gct.col_metadata_df)

        pd.util.testing.assert_frame_equal(expected_conn_gct.data_df, actual_conn_gct.data_df)
        pd.util.testing.assert_frame_equal(expected_conn_gct.row_metadata_df, actual_conn_gct.row_metadata_df)
        pd.util.testing.assert_frame_equal(expected_conn_gct.col_metadata_df, actual_conn_gct.col_metadata_df)

        os.remove(out_steep_path)
        os.remove(out_sip_path)


if __name__ == '__main__':
    unittest.main()
