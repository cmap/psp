import logging
import numpy as np
import os
import unittest

import broadinstitute_psp.utils.setup_logger as setup_logger
import cmapPy.pandasGEXpress.parse as parse
import broadinstitute_psp.external_query.external_query as eq

FUNCTIONAL_TESTS_DIR = "external_query/functional_tests"

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestExternalQuery(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        external_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_external_query_external.gct")
        cls.external_gct = parse.parse(external_gct_path)

        internal_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_external_query_internal.gct")
        cls.internal_gct = parse.parse(internal_gct_path)

        bg_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_external_query_bg.gct")
        cls.bg_gct = parse.parse(bg_gct_path)

    def test_replicates_in_external(self):
        # Test what happens when the external dataset has replicates
        (sim_gct, conn_gct) = eq.do_steep_and_sip(
            self.external_gct, self.internal_gct, self.bg_gct, "spearman",
            "ks_test", ["pert_id", "cell_id", "pert_time"],
            ["pert_id", "cell_id", "pert_time"])

        # Expected
        e_shape_conn_gct = (3, 3)
        e_dmso_v_etoposide = np.array([-0.4061, -0.5030, -0.3939, 0.9030, 0.8909, 0.9394, 0.5758, 0.4303, 0.5030])
        e_brdK08970894_v_brdK74148702 = -0.5556  # internal v. external
        e_brdA81177136_v_brdK87737963 = 0.3333
        e_conn_col_meta_columns = ["cell_id", "det_plate", "det_well",
                                   "pert_id", "pert_iname", "pert_time",
                                   "similarity_metric", "connectivity_metric"]
        e_col_meta_val1 = "A10:A11:A12"
        e_col_meta_val2 = "spearman"

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

        # Test some metadata
        self.assertItemsEqual(e_conn_col_meta_columns, conn_gct.col_metadata_df.columns)
        self.assertEqual(
            conn_gct.col_metadata_df.loc["BRD-K87737963:A375:3", "det_well"],
            e_col_meta_val1)
        self.assertEqual(
            conn_gct.col_metadata_df.loc["BRD-K87737963:A375:3", "similarity_metric"],
            e_col_meta_val2)

    def test_no_replicates_in_external(self):
        # Test what happens when the external dataset has no replicates
        (sim_gct, conn_gct) = eq.do_steep_and_sip(
            self.external_gct, self.internal_gct, self.bg_gct, "spearman",
            "ks_test", ["det_well"],
            ["pert_id", "cell_id", "pert_time"])

        # Expected
        e_shape_conn_gct = (9, 3)
        e_brdK08970894_v_a10 = -0.6667  # n.b. internal v. external
        e_brdA81177136_v_b5 = 0.3333
        e_conn_col_meta_columns = ["cell_id", "det_plate", "det_well",
                                   "pert_id", "pert_iname", "pert_time",
                                   "similarity_metric", "connectivity_metric"]

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

        # Test some metadata
        self.assertItemsEqual(e_conn_col_meta_columns, conn_gct.col_metadata_df.columns)
        self.assertEqual(
            conn_gct.col_metadata_df.loc["B3", "pert_iname"], "Etoposide")


if __name__ == '__main__':
    setup_logger.setup(verbose=True)
    unittest.main()
