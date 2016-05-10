import unittest
import numpy as np
import pandas as pd
# import robust_zscore

"""
FORMULA: x' = (x - row_median) / (row_mad * 1.4826)
"""

class TestRobustZscore(unittest.TestCase):
    def test_robust_zscore(self):
        my_df = pd.DataFrame([[0.3,0.2,0.8,0.7],
                             [0.1,-0.2,0.1,0.9],
                             [1.1,-0.3,np.nan,0.4]])
        e_df = pd.DataFrame([[0.3,0.2,0.8,0.7],
                             [0.1,-0.2,0.1,0.9],
                             [1.1,-0.3,np.nan,0.4]])

        row_medians = my_df.median(axis=1)
        row_mads = (my_df.sub(row_medians, axis=0)).abs().median(axis=1)

        # TODO(lev): continue here!

        print row_medians
        print row_mads


if __name__ == "__main__":
    unittest.main()
