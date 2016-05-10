import unittest
import logging
import utils.setup_logger as setup_logger
import numpy as np
import pandas as pd
import os
import qc_gct2pw

# N.B. The utils directory contains a directory called functional_tests
# that has the assets required below.

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Functional tests dir lives within the utils directory
FUNCTIONAL_TESTS_DIR = "utils/functional_tests"


class TestQCGct2Pw(unittest.TestCase):
    def test_assemble_output_df(self):
        plate_names = np.array(["plate1", "plate2", "plate3"])
        well_names = np.array(["A1", "A2", "A3"])
        e_col_names = ["plate_name", "well_name", "arg1", "arg2", "arg3"]

        out_df = qc_gct2pw.assemble_output_df(
            plate_names, well_names, arg3=[1,2,3], arg1=[7,8,9], arg2=[4,5,6])

        self.assertEqual(list(out_df.columns.values), e_col_names)

    def test_extract_plate_and_well_names(self):
        plate_field = "plate_name"
        well_field = "well_name"
        col_meta = pd.DataFrame([["a","b","c"],["a","e","f"]],
                                columns=["plate_name", "other", "well_name"])
        [out_plates, out_wells] = qc_gct2pw.extract_plate_and_well_names(
            col_meta, plate_field, well_field)

        self.assertTrue(np.array_equal(out_plates, ["a","a"]))
        self.assertTrue(np.array_equal(out_wells, ["c","f"]))

    def test_main(self):
        IN_GCT = os.path.join(FUNCTIONAL_TESTS_DIR, "qc_plate32.gct")
        OUT_PW = os.path.join(FUNCTIONAL_TESTS_DIR, "test_qc_gct2pw.pw")
        args_string = "{} {}".format(IN_GCT, OUT_PW)
        args = qc_gct2pw.build_parser().parse_args(args_string.split())

        qc_gct2pw.main(args)

        # Clean up
        os.remove(OUT_PW)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
