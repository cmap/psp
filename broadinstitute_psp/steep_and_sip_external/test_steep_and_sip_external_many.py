import glob
import logging
import os
import shutil
import unittest
import broadinstitute_psp.utils.setup_logger as setup_logger
import steep_and_sip_external_many as ssec

logger = logging.getLogger(setup_logger.LOGGER_NAME)

FUNCTIONAL_TESTS_DIR = "steep_and_sip_external/functional_tests"


class TestSteepAndSipExternalMany(unittest.TestCase):
    def test_read_config_file(self):

        (cells, internal_gct_dir, bg_gct_dir,
         fields_to_aggregate_for_external_profiles,
         fields_to_aggregate_for_internal_profiles,
         similarity_metric, connectivity_metric) = ssec.read_config_file("clue/psp_on_clue.cfg")

        self.assertEqual(cells[2], "MCF7")
        self.assertEqual(fields_to_aggregate_for_internal_profiles[1], "cell_id")
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

    # Slow (~30 sec.)
    @unittest.skipUnless(os.path.exists("/cmap/"), "/cmap/ needs to exist to run TestSteepAndSipExternalMany.test_main")
    def test_main(self):
        test_external = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_many_single_sample.gct")
        test_config = os.path.join(FUNCTIONAL_TESTS_DIR, "test_steep_and_sip_external_many.cfg")

        args_string = "-a {} -e {} -o {} -p {}".format(
            "GCP", test_external, FUNCTIONAL_TESTS_DIR, test_config)
        args = ssec.build_parser().parse_args(args_string.split())

        ssec.main(args)

        out_dirs = glob.glob(os.path.join(FUNCTIONAL_TESTS_DIR, "steep_and_sip_external_many_????_??_??*"))
        self.assertTrue(len(out_dirs) > 0)

        # Delete each output directory (could be more than 1 if previous cleanup failed)
        [shutil.rmtree(out_dir) for out_dir in out_dirs]

if __name__ == '__main__':
    setup_logger.setup(verbose=True)
    unittest.main()