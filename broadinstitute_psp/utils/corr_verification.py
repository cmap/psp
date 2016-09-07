# Script 2/2 to ensure that corr did the same thing between R and Python (June 2016)

import numpy as np
import pandas as pd

METHOD='spearman'
# my_df = pd.DataFrame([[3,2,np.nan], [5,9,3], [1,np.nan],[2,8,2], [4,1,8]])
my_df = pd.DataFrame.from_dict({"a":[3,5,1,2,4], "b":[2,9,np.nan,8,1], "c":[np.nan,3,5,2,8]})
print("in_df: \n{}".format(my_df))

# Compute correlation for whole df
out_cor = my_df.corr(method=METHOD)
print("all_corrs: \n{}".format(out_cor))

# Expect col1 v. col2 to be [3,5,2,4] v. [2,9,8,1]
col1_v_col2_cor = pd.Series([3,5,2,4]).corr(pd.Series([2,9,8,1]), method=METHOD)
print("a v. b: {}, {}".format(col1_v_col2_cor, np.isclose(col1_v_col2_cor, out_cor.iloc[0,1])))

# Expect col1 v. col3 to be [5,1,2,4] v. [3,5,2,8]
col1_v_col3_cor = pd.Series([5,1,2,4]).corr(pd.Series([3,5,2,8]), method=METHOD)
print("a v. c: {}, {}".format(col1_v_col3_cor, np.isclose(col1_v_col3_cor, out_cor.iloc[0,2])))

# Expect col2 v. col3 to be [9,8,1] v. [3,2,8]
col2_v_col3_cor = pd.Series([9,8,1]).corr(pd.Series([3,2,8]), method=METHOD)
print("b v. c: {}, {}".format(col2_v_col3_cor, np.isclose(col2_v_col3_cor, out_cor.iloc[1,2])))
