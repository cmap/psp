import glob
import logging
import os
import shutil
import unittest
import broadinstitute_psp.utils.setup_logger as setup_logger
import external_query_many as eqm

logger = logging.getLogger(setup_logger.LOGGER_NAME)

FUNCTIONAL_TESTS_DIR = "external_query/functional_tests"


class TestExternalQueryMany(unittest.TestCase):
    def test_read_config_file(self):

        (cells, internal_gct_dir, bg_gct_dir,
         fields_to_aggregate_for_internal_profiles,
         similarity_metric, connectivity_metric) = eqm.read_config_file("clue/psp_on_clue.cfg")

        self.assertEqual(cells[2], "MCF7")
        self.assertEqual(fields_to_aggregate_for_internal_profiles[1], "cell_id")
        self.assertEqual(connectivity_metric, "ks_test")

    def test_write_failure(self):
        out_file = os.path.join(FUNCTIONAL_TESTS_DIR, "failure.txt")

        try:
            tuple()[0]
        except Exception:
            eqm.write_failure(out_file, "Started way back when!")

        self.assertTrue(os.path.exists(out_file))
        os.remove(out_file)

    def test_write_success(self):
        out_file = os.path.join(FUNCTIONAL_TESTS_DIR, "success.txt")

        eqm.write_success(out_file, "Started now!")

        self.assertTrue(os.path.exists(out_file))
        os.remove(out_file)

    # Slow (~30 sec.)
    @unittest.skipUnless(os.path.exists("/cmap/"), "/cmap/ needs to exist to run TestExternalQueryMany.test_main")
    def test_main(self):
        test_external = os.path.join(FUNCTIONAL_TESTS_DIR, "test_external_query_many_single_sample.gct")
        test_config = os.path.join(FUNCTIONAL_TESTS_DIR, "test_external_query_many.cfg")

        args_string = "-a {} -e {} -o {} -p {} -fae {} -i".format(
            "GCP", test_external, FUNCTIONAL_TESTS_DIR, test_config, "pert_id cell_id")
        args = eqm.build_parser().parse_args(args_string.split())
        logger.info(args)

        this_uuid = eqm.main(args)

        out_dir = os.path.join(FUNCTIONAL_TESTS_DIR, "external_query_many_" + this_uuid)

        self.assertTrue(os.path.exists(os.path.join(out_dir, "CONCATED_CONN.gct")))
        self.assertTrue(os.path.exists(os.path.join(out_dir, "INTROSPECT_CONN.gct")))

        # Clean up
        shutil.rmtree(out_dir)


if __name__ == '__main__':
    setup_logger.setup(verbose=True)
    unittest.main()