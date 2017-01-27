import logging
import pandas as pd
import unittest
import network_query_of_gct_igraph as nq
import broadinstitute_psp.utils.setup_logger as setup_logger

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestNetworkQueryOfGctIgraph(unittest.TestCase):

    def test_make_color_dict(self):
        meta_df = pd.DataFrame(
            [["A375", "A"],
             ["A375", "B"],
             ["A375", "A"]],
            index=["a", "b", "c"],
            columns=["cell", "moa"])
        e_color_dict = {'A': (1, 0, 1, 1), 'B': (0, 1, 1, 1)}
        color_dict = nq.make_color_dict(meta_df, "moa")
        self.assertEqual(color_dict, e_color_dict)


if __name__ == '__main__':
    setup_logger.setup(verbose=True)
    logger.info("here")
    unittest.main()
