import unittest
import logging
import utils.setup_logger as setup_logger
import numpy as np
import gct2pw

"""
This code should be run from broadinstitute.psp/utils.

The utils directory contains a directory called functional_tests
that has the assets required below.

"""

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestGct2Pw(unittest.TestCase):
    def test_assemble_output_df(self):
        plate_names = np.array(["plate1", "plate2", "plate3"])
        well_names = np.array(["A1", "A2", "A3"])
        e_col_names = ["plate_name", "well_name", "arg1", "arg2", "arg3"]

        out_df = gct2pw.assemble_output_df(
            plate_names, well_names, arg3=[1,2,3], arg1=[7,8,9], arg2=[4,5,6])

        self.assertEqual(list(out_df.columns.values), e_col_names)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
