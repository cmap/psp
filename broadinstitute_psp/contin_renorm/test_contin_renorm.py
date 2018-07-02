import unittest
import logging
import os
import numpy as np
import pandas as pd

import cmapPy.pandasGEXpress as GCToo
import cmapPy.pandasGEXpress.parse as parse
import broadinstitute_psp.utils.setup_logger as setup_logger
import contin_renorm as renorm

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)

# Use functional tests assets from the tear directory
FUNCTIONAL_TESTS_DIR = os.path.join("functional_tests")


class test_contin_renorm(unittest.TestCase):
    def test_main(self):
        in_gct_path = os.path.join(FUNCTIONAL_TESTS_DIR, "test_renorm_main.gct")
        out_name = os.path.join(FUNCTIONAL_TESTS_DIR, "test_renorm_out.gct")

        # update the args string
        args_string = ("-i {} -o {}").format(in_gct_path, out_name)
        args = renorm.build_parser().parse_args(args_string.split())
        renorm.contin_renorm_main(args)

        # Read in result
        out_gct = parse.parse(out_name)

        e_values = np.array(
            [[-0.41, -0.13, 0.07, 0.09, 0.18, 0.24, 0.08],
             [0.40, 0.11, 0.06, -0.11, -0.22, -0.26, -0.09],
             [0.40, -0.40, 0.30, -0.20, 0.05, -0.10, 0.10],
             [0.10, 0.06, -0.07, 0.05, -0.09, 0.08, 0.10]])
        self.assertTrue(np.allclose(e_values, out_gct.data_df, atol=1e-2))

        # Clean up
        os.remove(out_name)
        
        
    def test_to_log_renorm(self):
        
        slopes = pd.Series([-0.3, 0.1, -0.1, 0.3])
        is_log_renormed_return = renorm.is_log_renormed(slopes, 0.2)
        is_log_renormed_expected = pd.Series([True, False, False, True])
        
        self.assertTrue((is_log_renormed_return == is_log_renormed_expected).all())
        
    
    def test_tot_samp_offsets(self):
        
        df_in = pd.DataFrame([[1,     -1,   0, 1],
                              [2,      0,  -2, 2],
                              [0,      0,   0, 0],
                              [-0.5, 0.5, 0.5, 0]])
        return_tot_samp_offsets = renorm.calc_tot_samp_offsets(df_in)
        expect_tot_samp_offsets = pd.Series([3.5, 1.5, 2.5, 3])
        
        self.assertTrue(np.allclose(return_tot_samp_offsets,
                                    expect_tot_samp_offsets,
                                    atol=1e-6))
        
    
    def test_calc_out_mat(self):
        
        df_in = pd.DataFrame([[1, 2, 3, 4],
                              [4, 3, 2, 1],
                              [0, 0, 0, 0],
                              [1, 1, 1, 1]])
        offset_in = pd.DataFrame([[0, 0, 0, 0],
                                  [3, 2, 1, 0],
                                  [0, 1, 2, 3],
                                  [0, 0, 0, 0]])
        return_out_df = renorm.calc_out_mat(df_in, offset_in)
        
        expect_out_df = pd.DataFrame([[1, 2, 3, 4],
                                      [1, 1, 1, 1],
                                      [0, -1, -2, -3],
                                      [1, 1, 1, 1]])
        
        self.assertTrue(np.allclose(return_out_df, expect_out_df,
                                    atol=1e-6))
            
    
    def test_calc_pep_samp_offsets(self):
        
        data_df_in = pd.DataFrame([[0.8, 0.6, 0.5, 0.36],
                                   [1, 1, 1, 1]])
        row_metadata_df_in = pd.DataFrame([[True, True],
                                           [False, False]],
                                          columns = ["is_log_renormed", "whatever"])
        es_in = pd.Series([0.2, 0.5, 0.6, 1.0])
        pep_y_offsets_in = pd.Series([0.4, 0])
        fit_params_in = pd.DataFrame([[(1, 1), (1, 1)],
                                      [(1, 1), (1, 1)]],
                                     columns = ["deg1", "log"])

        func_return = renorm.calc_pep_samp_offsets(data_df_in, row_metadata_df_in,
                                                   es_in, fit_params_in,
                                                   pep_y_offsets_in)
        
        expected_return = pd.DataFrame([[0.85, 0.78, 0.75, 0.67],
                                        [0,    0,    0,    0]])
        
        self.assertTrue(np.allclose(expected_return, func_return, atol=1e-2))
        
            
    def test_calc_fit(self):
        
        data_df_in = pd.DataFrame([[0.8, 0.6, 0.5, 0.36],
                                   [0.9, 1, 1, 1]])
        es_in = pd.Series([0.2, 0.5, 0.6, 1.0])
        pep_y_offsets_in = pd.Series([0.1, 0.1])
        
        func_return = renorm.calc_fit(data_df_in, es_in, pep_y_offsets_in)

        expect_return = pd.DataFrame([[[-0.54, 0.88], (1.6, 1.66)],
                                      [[0.11, 0.91],  (1.8, 0.03)]],
                                     columns=["deg1", "log"])
        
        for row_idx, vals in expect_return.iterrows():
            for col_idx, vals in expect_return.iteritems():
                self.assertTrue(np.allclose(expect_return.loc[row_idx, col_idx],
                                            func_return.loc[row_idx, col_idx],
                                            atol=1e-2))
        
        
    def test_make_y(self):
        
        x_in = pd.Series([0.1, 0.3, 0.5, 0.8])
        deg_model = [1, 1]
        log_model = (1, 1)
        
        deg_return = renorm.make_y(x_in, deg_model)
        log_return = renorm.make_y(x_in, log_model, 1)
        
        expect_deg_return = pd.Series([1, 1, 1, 1])
        expect_log_return = pd.Series([1.47, 1.42, 1.37, 1.31])
        
        self.assertTrue(np.allclose(deg_return, expect_deg_return, atol=1e-2))
        self.assertTrue(np.allclose(log_return, expect_log_return, atol=1e-2))
        
        
    def test_calc_y_offsets(self):
        
        df_in = pd.DataFrame([[ 1,  2,  3,  4,  5],
                              [ 4,  3,  2,  1,  0],
                              [ 0,  0,  0,  0,  0],
                              [-1, -1, -1, -1, -1]])
        es = pd.Series([1, 0.2, 0.3, 0.4, 0.5])
        
        return_y_offsets = renorm.calc_y_offsets(df_in, es)
        expect_y_offsets = pd.Series([1, 4, 0, -1])
        
        self.assertTrue(np.allclose(return_y_offsets, expect_y_offsets,
                                            atol=1e-6))
        
        return_y_offsets = renorm.calc_y_offsets(df_in, es, top_frac=1.0)
        expect_y_offsets = pd.Series([3, 2, 0, -1])
        
        self.assertTrue(np.allclose(return_y_offsets, expect_y_offsets,
                                            atol=1e-6))
    
    

if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()
