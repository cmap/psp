import unittest
import steep_and_sip_external_many as ssec

# HOWTO: run from broadinstitute_psp directory


class TestSteepAndSipExternalMany(unittest.TestCase):
    def test_read_config_file(self):

        (cells, internal_gct_dir, bg_gct_dir,
         fields_to_aggregate_for_external_profiles,
         fields_to_aggregate_for_internal_profiles,
         similarity_metric, connectivity_metric) = ssec.read_config_file("clue/psp_on_clue.cfg")

        self.assertEqual(cells[2], "MCF7")
        self.assertEqual(fields_to_aggregate_for_external_profiles[1], "cell_id")
        self.assertEqual(connectivity_metric, "ks_test")


if __name__ == '__main__':
    unittest.main()
