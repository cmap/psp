import logging
import numpy as np
import os
import unittest

import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_cmap.io.pandasGEXpress.parse_gct as pg
import broadinstitute_psp.steep_and_sip_external.steep_and_sip_external as sse

FUNCTIONAL_TESTS_DIR = "steep_and_sip_external/functional_tests"

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestSteepAndSipExternal(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        external_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_external.gct")
        cls.external_gct = pg.parse(external_gct_path, convert_neg_666=False, make_multiindex=True)

        internal_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_internal.gct")
        cls.internal_gct = pg.parse(internal_gct_path, convert_neg_666=False, make_multiindex=True)

        bg_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_bg.gct")
        cls.bg_gct = pg.parse(bg_gct_path, convert_neg_666=False, make_multiindex=True)

    def test_replicates_in_external(self):
        # Test what happens when the external dataset has replicates
        (sim_gct, conn_gct) = sse.do_steep_and_sip(
            self.external_gct, self.internal_gct, self.bg_gct, "spearman",
            "ks_test", ["pert_id", "cell_id", "pert_time"],
            ["pert_id", "cell_id", "pert_time"])

        # Expected
        e_shape_conn_gct = (3, 3)
        e_dmso_v_etoposide = np.array([-0.4061, -0.5030, -0.3939, 0.9030, 0.8909, 0.9394, 0.5758, 0.4303, 0.5030])
        e_brdK08970894_v_brdK74148702 = -0.5556  # internal v. external
        e_brdA81177136_v_brdK87737963 = 0.3333
        e_conn_col_meta_columns = ["pert_id", "cell_id", "pert_time", "connectivity_metric"]

        # Test that shape correct
        self.assertItemsEqual(e_shape_conn_gct, conn_gct.data_df.shape)

        # Test a few similarities
        np.testing.assert_allclose(
            e_dmso_v_etoposide, sim_gct.data_df.loc[
                sim_gct.row_metadata_df.pert_iname == "DMSO",
                sim_gct.col_metadata_df.pert_iname == "Etoposide"].values.flatten(), atol=1e-4)

        # Test a few connectivities
        self.assertAlmostEqual(e_brdK08970894_v_brdK74148702, float(
            conn_gct.data_df.loc[conn_gct.row_metadata_df.pert_id == "BRD-K08970894",
                                 conn_gct.col_metadata_df.pert_id == "BRD-K74148702"].values), places=4)
        self.assertAlmostEqual(e_brdA81177136_v_brdK87737963, float(
            conn_gct.data_df.loc[conn_gct.row_metadata_df.pert_id == "BRD-A81177136",
                                 conn_gct.col_metadata_df.pert_id == "BRD-K87737963"].values), places=4)
        self.assertItemsEqual(e_conn_col_meta_columns, conn_gct.col_metadata_df.columns)

    def test_no_replicates_in_external(self):
        # Test what happens when the external dataset has no replicates
        (sim_gct, conn_gct) = sse.do_steep_and_sip(
            self.external_gct, self.internal_gct, self.bg_gct, "spearman",
            "ks_test", ["det_well"],
            ["pert_id", "cell_id", "pert_time"])

        # Expected
        e_shape_conn_gct = (9, 3)
        e_brdK08970894_v_a10 = -0.6667  # n.b. internal v. external
        e_brdA81177136_v_b5 = 0.3333
        e_conn_col_meta_columns = ["det_well", "connectivity_metric"]

        # Test that shape correct
        self.assertItemsEqual(e_shape_conn_gct, conn_gct.data_df.shape)

        # Test a few connectivities
        self.assertAlmostEqual(
            e_brdK08970894_v_a10, float(conn_gct.data_df.loc[
                conn_gct.row_metadata_df.pert_id == "BRD-K08970894",
                conn_gct.col_metadata_df.det_well == "A10"].values.flatten()), places=4)
        self.assertAlmostEqual(
            e_brdA81177136_v_b5, float(conn_gct.data_df.loc[
                conn_gct.row_metadata_df.pert_id == "BRD-A81177136",
                conn_gct.col_metadata_df.det_well == "B5"].values.flatten()), places=4)
        self.assertItemsEqual(e_conn_col_meta_columns, conn_gct.col_metadata_df.columns)


if __name__ == '__main__':
    setup_logger.setup(verbose=True)
    unittest.main()
