import unittest
import logging
import setup_logger as setup_logger
import dry

logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestDry(unittest.TestCase):
    def test_update_prov_code(self):
#         (new_entry, existing_prov_code):
        pass


	def test_filter_samples(self):
# 	    filter_samples(df, sample_pct_cutoff)
	    # Return df (potentially of different size than original df)

    def test_filter_probes(self):
#         filter_probes(df, probe_pct_cutoff, probe_sd_cutoff)
        # Identify rows manually labeled for rejection
        # Return df (potentially of different size than original df)
        pass

	def test_main(self):
	    dry.main()


if __name__ == "__main__":
	setup_logger.setup(verbose=True)

	unittest.main()