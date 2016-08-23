import unittest
import logging
import utils.setup_logger as setup_logger
import os
import numpy as np
import pandas as pd
import scipy.stats as stats

import parse_gctoo
import GCToo
import old_steep

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Functional tests dir lives within the dry directory
FUNCTIONAL_TESTS_DIR = "steep/functional_tests"

## TODO(lev): does not currently work

class TestSteep(unittest.TestCase):

    def test_compute_similarity_against_R_code(self):
        INPUT_GCT_PATH = os.path.join(FUNCTIONAL_TESTS_DIR, "PC3_QCNORM_n189x59.gct")
        JJ_PC3_SIM_PATH = os.path.join(FUNCTIONAL_TESTS_DIR, "R_PC3_SIM_n189x189.csv")

        # Read intest_steep.py PC3 data and JJ's R version of the similarity matrix
        pc3_qcnorm_gct = parse_gctoo.parse(INPUT_GCT_PATH)
        jj_pc3_sim_df = pd.read_csv(JJ_PC3_SIM_PATH, index_col=0)

        # Compute similarity
        pc3_sim_df = old_steep.compute_similarity_matrix(pc3_qcnorm_gct.data_df, method="spearman")

        self.assertTrue(np.allclose(pc3_sim_df, jj_pc3_sim_df), (
            "\npc3_sim_df:\n{}\njj_pc3_sim_df:\n{}".format(pc3_sim_df, jj_pc3_sim_df)))



if __name__ == "__main__":
    setup_logger.setup(verbose=False)
    unittest.main()