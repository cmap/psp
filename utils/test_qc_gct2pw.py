import unittest
import logging
import utils.setup_logger as setup_logger
import numpy as np
import pandas as pd
import os
import qc_gct2pw

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestQCGct2Pw(unittest.TestCase):

    def test_assemble_output_df(self):
        plate_names = np.array(["plate1", "plate2", "plate3"])
        well_names = np.array(["A1", "A2", "A3"])
        metadata_dict = {"arg3": [1, 2, 3], "arg1": [7, 8, 9], "arg2": [4, 5, 6]}
        e_col_names = ["plate_name", "well_name", "arg1", "arg2", "arg3"]

        out_df = qc_gct2pw.assemble_output_df(
            plate_names, well_names, metadata_dict)
        self.assertEqual(list(out_df.columns.values), e_col_names)

    def test_extract_plate_and_well_names(self):
        plate_field = "plate_name"
        well_field = "well_name"
        col_meta = pd.DataFrame(
            [["a", "b", "c"], ["a", "e", "f"]],
            columns=["plate_name", "other", "well_name"])

        [out_plates, out_wells] = qc_gct2pw.extract_plate_and_well_names(
            col_meta, plate_field, well_field)
        self.assertTrue(np.array_equal(out_plates, ["a", "a"]))
        self.assertTrue(np.array_equal(out_wells, ["c", "f"]))

    def test_undo_log_transform_if_needed(self):
        data_df = pd.DataFrame([[1,2,4],[2,3,4]])
        prov_code = ["L2X"]
        e_out_df = pd.DataFrame([[2,4,16],[4,8,16]])

        out_df = qc_gct2pw.undo_log_transform_if_needed(data_df, prov_code)

        self.assertTrue(np.allclose(out_df, e_out_df),
                        "out_df is incorrect: {}".format(out_df))

    def test_main(self):
        in_gct = "utils/functional_tests/test_qc.gct"
        out_pw = "utils/functional_tests/test_qc_gct2pw_output.pw"
        args_string = "{} {}".format(in_gct, out_pw)
        args = qc_gct2pw.build_parser().parse_args(args_string.split())

        qc_gct2pw.main(args)

        # Verify that pw file exists
        self.assertTrue(os.path.exists(out_pw))

        # Remove pw file
        os.remove(out_pw)

if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
