import os
import unittest
import steep_and_sip_external_many as ssec

FUNCTIONAL_TESTS_DIR = "steep_and_sip_external/functional_tests"


class TestSteepAndSipExternalMany(unittest.TestCase):
    def test_read_config_file(self):

        (cells, internal_gct_dir, bg_gct_dir,
         fields_to_aggregate_for_external_profiles,
         fields_to_aggregate_for_internal_profiles,
         similarity_metric, connectivity_metric) = ssec.read_config_file("clue/psp_on_clue.cfg")

        self.assertEqual(cells[2], "MCF7")
        self.assertEqual(fields_to_aggregate_for_external_profiles[1], "cell_id")
        self.assertEqual(connectivity_metric, "ks_test")

    def test_write_failure(self):
        out_file = os.path.join(FUNCTIONAL_TESTS_DIR, "failure.txt")

        try:
            tuple()[0]
        except Exception:
            ssec.write_failure(out_file)

        self.assertTrue(os.path.exists(out_file))
        os.remove(out_file)

    def test_write_success(self):
        out_file = os.path.join(FUNCTIONAL_TESTS_DIR, "success.txt")

        ssec.write_success(out_file)

        self.assertTrue(os.path.exists(out_file))
        os.remove(out_file)

    def test_main(self):
        test_external = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_many_single_sample.gct")

        args_string = "-a {} -e {} -o {}".format("GCP", test_external, FUNCTIONAL_TESTS_DIR)
        args = ssec.build_parser().parse_args(args_string.split())

        ssec.main(args)

if __name__ == '__main__':
    unittest.main()
