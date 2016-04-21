import unittest
import logging
import utils.setup_logger as setup_logger
import qc_gct_to_pw

"""
This code should be run from broadinstitute.psp/utils.

The utils directory contains a directory called functional_tests
that has the assets required below.

"""

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestQCGctToPw(unittest.TestCase):
    def test_main(self):
        gct_path = ""

        qc_gct_to_pw.main()


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
