import unittest
import logging
import os
import tempfile
import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.utils.psp_utils as psp_utils
import broadinstitute_psp.utils.config_converter as config_converter

logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestConfigConverter(unittest.TestCase):

    def test_convert_config_to_json(self):
        # Differential parameters, returns JSON
        json = config_converter.convert_config_to_json("gcp", "functional_tests/test_config_converter_diff.cfg")
        expected_json = "{'gcp_sample_frac_cutoff':0.7,'gcp_probe_frac_cutoff':0.7}"
        self.assertEqual(json, expected_json)

        # No differential parameters, returns None
        json = config_converter.convert_config_to_json("gcp", "functional_tests/test_config_converter_no_diff.cfg")
        self.assertIsNone(json)

    def test_convert_json_to_config(self):
        config_dir = tempfile.mkdtemp()
        save_file_path = os.path.join(config_dir, "test.cfg")
        assay = "gcp"
        # probe_frac_cutoff == default value
        json = '""{""samplePctCutoff"":0.7,""probePctCutoff"":0.5}""'
        config_converter.convert_json_to_config(assay, json, save_file_path)

        (_, _, custom_parameters_from_read) = psp_utils.read_config_file(save_file_path)
        differential_params = config_converter.check_custom_parameters_against_defaults(assay, custom_parameters_from_read)

        expected_params = {"gcp_sample_frac_cutoff":"0.7"}
        self.assertEqual(differential_params, expected_params)

    def test_convert_gct_to_config(self):
        # Differential params
        config_dir = tempfile.mkdtemp()
        save_file_path = os.path.join(config_dir, "test.cfg")

        gct_path = "functional_tests/test_pr_processing_params.gct"
        assay = "p100"
        diff_params = config_converter.convert_gct_to_config(assay, gct_path, save_file_path)
        self.assertIsNotNone(diff_params)

        (_, _, custom_parameters_from_read) = psp_utils.read_config_file(save_file_path)
        differential_params = config_converter.check_custom_parameters_against_defaults(assay, custom_parameters_from_read)

        expected_params = {"p100_dist_sd_cutoff": "6"}
        self.assertEqual(differential_params, expected_params)

        # No differential params, pr_processing_params empty {}
        gct_path = "functional_tests/test_pr_processing_params_empty.gct"
        diff_params = config_converter.convert_gct_to_config(assay, gct_path, save_file_path)
        self.assertIsNone(diff_params)

        # No differential params, pr_processing_params column DNE
        gct_path = "functional_tests/test_no_pr_processing_params.gct"
        diff_params = config_converter.convert_gct_to_config(assay, gct_path, save_file_path)
        self.assertIsNone(diff_params)

        # No differential params, pr_processing_params contains base_histone_normalization_mapping only
        gct_path = "functional_tests/test_pr_processing_params_base_histone_normalization_mapping.gct"
        diff_params = config_converter.convert_gct_to_config(assay, gct_path, save_file_path)
        self.assertIsNone(diff_params)


    def test_create_dict_from_pseudojson(self):
        pseudoJSON = '""{""samplePctCutoff"":0.7,""probePctCutoff"":0.5}""'
        dictionary_returned = config_converter.create_dict_from_pseudojson(pseudoJSON)
        dictionary_expected = {"samplePctCutoff":"0.7",
                               "probePctCutoff":"0.5"}
        self.assertEqual(dictionary_expected, dictionary_returned)

    def test_populate_background_parameters(self):
        # single custom params with differential value assigned
        custom_params_unmapped = {"probe_sd_cutoff": "5"}
        assay = "p100"
        mapped_params = config_converter.map_R_params(assay, custom_params_unmapped)
        returned_full_dict = config_converter.populate_background_parameters(mapped_params)
        full_dict= {
            'p100_probe_sd_cutoff': '5',
            'offset_bounds': '(-7, 7)',
            'p100_probe_frac_cutoff': '0.9',
            'gcp_sample_frac_cutoff': '0.5',
            'p100_sample_frac_cutoff': '0.8',
            'gcp_probe_sd_cutoff': '4',
            'gcp_probe_frac_cutoff': '0.5'
        }
        self.assertEqual(returned_full_dict, full_dict)

    def test_check_custom_parameters_against_defaults(self):
        # No differential params from pseudoJSON
        assay = "p100"
        pseudoJSON = '""{""samplePctCutoff"":0.8,""probePctCutoff"":0.9}""'
        no_diff_dict = config_converter.create_dict_from_pseudojson(pseudoJSON)
        differential_parameters = config_converter.check_custom_parameters_against_defaults(assay,no_diff_dict, json=True)
        self.assertIsNone(differential_parameters)

        # Differential params from pseudoJSON
        pseudoJSON = '""{""samplePctCutoff"":0.7,""probePctCutoff"":0.5}""'
        diff_dict = config_converter.create_dict_from_pseudojson(pseudoJSON)
        differential_parameters = config_converter.check_custom_parameters_against_defaults(assay, diff_dict, json=True)
        expected_params = {"p100_sample_frac_cutoff":"0.7", "p100_probe_frac_cutoff":"0.5"}
        self.assertEqual(differential_parameters, expected_params)

        # No differential params, json=False
        no_diff_dict = {
            'p100_probe_sd_cutoff': '3',
            'offset_bounds': '(-7, 7)',
            'p100_probe_frac_cutoff': '0.9',
            'gcp_sample_frac_cutoff': '0.5',
            'p100_sample_frac_cutoff': '0.8',
            'gcp_probe_sd_cutoff': '4',
            'gcp_probe_frac_cutoff': '0.5'
        }
        assay = "p100"
        differential_parameters = config_converter.check_custom_parameters_against_defaults(assay, no_diff_dict)
        self.assertIsNone(differential_parameters)

        # Differential params, json=False
        diff_dict = {
            'p100_probe_sd_cutoff': '4',
            'offset_bounds': '(-7, 7)',
            'p100_probe_frac_cutoff': '0.9',
            'gcp_sample_frac_cutoff': '0.5',
            'p100_sample_frac_cutoff': '0.8',
            'gcp_probe_sd_cutoff': '4',
            'gcp_probe_frac_cutoff': '0.5'
        }
        assay = "p100"
        differential_parameters = config_converter.check_custom_parameters_against_defaults(assay, diff_dict)
        expected_params = {"p100_probe_sd_cutoff":"4"}
        self.assertEqual(differential_parameters, expected_params)

    def test_map_R_params(self):
        # Match to R param name
        dict = {"samplePctCutoff":0.7}
        assay= "p100"
        returned_params = config_converter.map_R_params(assay, dict)
        expected_params = {"p100_sample_frac_cutoff":0.7}
        self.assertEqual(returned_params,expected_params)

        # Match psp param name without assay pre-pended
        dict = {"sample_frac_cutoff":0.7}
        assay = "p100"
        returned_params = config_converter.map_R_params(assay, dict)
        expected_params = {"p100_sample_frac_cutoff":0.7}
        self.assertEqual(returned_params, expected_params)

        # No match, valid psp param check passes
        dict = {"p100_sample_frac_cutoff":0.7}
        assay = "p100"
        returned_params = config_converter.map_R_params(assay, dict)
        expected_params = dict
        self.assertEqual(returned_params, expected_params)

        # No match, valid psp param check fails

        dict = {"p100_probe_frac_cutoff":0.8,"sample_frac_cutoff_p100":0.7}
        assay = "p100"
        returned_params = config_converter.map_R_params(assay, dict)
        self.assertIsNone(returned_params)

    def test_write_config(self):
        config_dir = tempfile.mkdtemp()
        save_file_path = os.path.join(config_dir, "test.cfg")
        custom_params = {"p100_sample_frac_cutoff":"0.1"}
        config_converter.write_config(custom_params, save_file_path)

        (_,_, custom_parameters_from_read) = psp_utils.read_config_file(save_file_path)
        differential_params = config_converter.check_custom_parameters_against_defaults("p100", custom_parameters_from_read)
        self.assertEqual(differential_params, custom_params)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()