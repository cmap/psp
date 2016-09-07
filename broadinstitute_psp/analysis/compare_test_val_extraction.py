import pandas as pd
import numpy as np

df_index = pd.MultiIndex.from_arrays(
    [["A", "B", "A", "B", "C", "C"], [1, 2, 3, 4, 5, 6]], names=["group", "id"])
df = pd.DataFrame(
    [[1.0, 0.5, 1.0, -0.4, 1.1, -0.6],
     [0.5, 1.0, 1.2, -0.8, -0.9, 0.4],
     [1.0, 1.2, 1.0, 0.1, 0.3, 1.3],
     [-0.4, -0.8, 0.1, 1.0, 0.5, -0.2],
     [1.1, -0.9, 0.3, 0.5, 1.0, 0.7],
     [-0.6, 0.4, 1.3, -0.2, 0.7, 1.0]],
    index=df_index, columns=df_index)


def extract_test_vals2(query, target, multi_index_level_name, df):
    # Get indices for what you want to extract and then extract all at once
    idxs = [[i, j] for i in range(len(df)) for j in range(len(df)) if i < j and (
        df.index.get_level_values(multi_index_level_name)[i] == target and (
            df.columns.get_level_values(multi_index_level_name)[j] == query))]

    return df.values[tuple(np.transpose(idxs))]


def extract_test_vals(query, target, multi_index_level_name, df):
    # Extract entries where target is in the rows and query is in the columns
    target_in_rows_query_in_cols_df = df.loc[
        df.index.get_level_values(multi_index_level_name) == target,
        df.columns.get_level_values(multi_index_level_name) == query]
    mask = np.triu(np.ones(target_in_rows_query_in_cols_df.shape), k=1).astype(np.bool)
    vals_with_nans = target_in_rows_query_in_cols_df.where(mask).values.flatten()
    vals = vals_with_nans[~np.isnan(vals_with_nans)]

    return vals

# Expected values
e_A_B_vals = np.sort([0.5, -0.4, 1.2, 0.1])
e_A_C_vals = np.sort([1.1, 0.3, -0.6, 1.3])
e_C_A_vals = np.sort([1.1, 0.3, -0.6, 1.3])
e_A_A_vals = np.sort([1.0])
e_C_C_vals = np.sort([0.7])

# # Sort because order doesn't matter
# assert np.allclose(np.sort(extract_test_vals("A", "B", "group", df)), e_A_B_vals)
# assert np.allclose(np.sort(extract_test_vals("A", "C", "group", df)), e_A_C_vals)
# assert np.allclose(np.sort(extract_test_vals("C", "A", "group", df)), e_C_A_vals)
# assert np.allclose(np.sort(extract_test_vals("A", "A", "group", df)), e_A_A_vals)
# assert np.allclose(np.sort(extract_test_vals("C", "C", "group", df)), e_C_C_vals)

print extract_test_vals2("A", "B", "group", df)

assert np.allclose(np.sort(extract_test_vals2("A", "B", "group", df)), e_A_B_vals)
assert np.allclose(np.sort(extract_test_vals2("A", "C", "group", df)), e_A_C_vals)
assert np.allclose(np.sort(extract_test_vals2("C", "A", "group", df)), e_C_A_vals)
assert np.allclose(np.sort(extract_test_vals2("A", "A", "group", df)), e_A_A_vals)
assert np.allclose(np.sort(extract_test_vals2("C", "C", "group", df)), e_C_C_vals)

## Test speed
import time

# Method 1
start1 = time.time()
for ii in range(1000):
    out = extract_test_vals("A", "B", "group", df)
elapsed1 = time.time() - start1
print elapsed1

# Method 2
start2 = time.time()
for ii in range(1000):
    out2 = extract_test_vals2("A", "B", "group", df)
elapsed2 = time.time() - start2
print elapsed2