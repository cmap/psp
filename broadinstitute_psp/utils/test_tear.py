import logging
import numpy as np
import os
import pandas as pd
import unittest

import broadinstitute_psp.utils.setup_logger as setup_logger
import tear

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Use functional tests assets from the dry directory
FUNCTIONAL_TESTS_DIR = "dry/functional_tests"


class TestTear(unittest.TestCase):
    def test_robust_zscore(self):
        my_df = pd.DataFrame([[0.3, 0.2, 0.8, 0.7],
                             [0.1, -0.2, 0.1, 0.9],
                             [1.1, -0.3, np.nan, 0.4]])
        e_df = pd.DataFrame([[-0.54, -0.81, 0.81, 0.54],
                             [0.0, -1.35, 0.0, 3.60],
                             [0.67, -0.67, np.nan, 0.0]])

        out_df = tear.robust_zscore(my_df)
        self.assertTrue(np.allclose(out_df, e_df, atol=1e-2, equal_nan=True))

    # TODO(lev): fix up
    # def test_main(self):
    #     INPUT_GCT_PATH = os.path.join(FUNCTIONAL_TESTS_DIR, "p100_prm_plate29_3H.gct")
    #     OUT_NAME = "test_tear_p100_output.gct"
    #
    #     args_string = ("{} {} -out_name {}").format(
    #         INPUT_GCT_PATH, FUNCTIONAL_TESTS_DIR, OUT_NAME)
    #     args = tear.build_parser().parse_args(args_string.split())
    #     out_gct = tear.main(args)
    #
    #     # print
    #
    #     self.assertAlmostEqual(out_gct.data_df.iloc[0,0], -1.20, places=2)
    #     self.assertEqual(out_gct.col_metadata_df["provenance_code"].iloc[0],
    #                      "PRM+L2X+ZSC")
    #
    #     # Clean up
    #     os.remove(os.path.join(FUNCTIONAL_TESTS_DIR, OUT_NAME))


if __name__ == "__main__":
    setup_logger.setup()
    unittest.main()
