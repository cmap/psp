import unittest
import logging
import os
import pandas as pd

import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.utils.psp_utils as utils

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestPspUtils(unittest.TestCase):

    def test_read_gct_and_config_file(self):
        config_path = "utils/functional_tests/test.cfg"
        assert os.path.exists(config_path), (
            "Config file cannot be found: {}".format(config_path))

        in_gct = "utils/functional_tests/test_p100.gct"
        assert os.path.exists(in_gct), (
            "in_gct cannot be found. in_gct: {}".format(in_gct))

        e_data_df_shape = (96, 96)
        e_nan_values = '["NA", "N/A"]'
        e_gcp_assays = '["GCP", "GR1"]'

        # Happy path
        (out_gct, config_io, config_meta, config_params) = (
            utils.read_gct_and_config_file(in_gct, config_path))

        self.assertEqual(out_gct.data_df.shape, e_data_df_shape, (
            "out_gct.data_df.shape is incorrect: {}").format(out_gct.data_df.shape))
        self.assertEqual(config_io["nan_values"], e_nan_values, (
            "config_io['nan_values'] is incorrect: {}".format(config_io["nan_values"])))
        self.assertEqual(config_meta["gcp_assays"], e_gcp_assays, (
            "config_meta['gcp_assays'] is incorrect: {}".format(config_meta["gcp_assays"])))
        self.assertEqual(config_params, {}, (
            "config_params should be empty: {}".format(config_params)))

    def test_extract_prov_code(self):
        col_meta_df = pd.DataFrame.from_dict(
            {"foo": ["a", "b", "c"],
             "prov_field": ["PRM+L2X", "PRM+L2X", "PRM+L2X"]})
        e_prov_code = ["PRM", "L2X"]

        prov_code = utils.extract_prov_code(col_meta_df, "prov_field", "+")
        self.assertEqual(e_prov_code, prov_code, (
            "prov_code is incorrect: {}").format(prov_code))

if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
