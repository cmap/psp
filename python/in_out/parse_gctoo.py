"""Reads in a gct file as a gctoo object.

The main method is parse. parse_into_3_dfs creates the row
metadata, column metadata, and data dataframes, while the
assemble_multi_index_df method in GCToo.py assembles them.
N.B. Only supports v1.3 gct files.

Example GCT:
#1.3
96  36  9 15

96 = number of data rows
36 = number of data columns
9 = number of row metadata fields (+1 for the 'id' column -- first column)
15 = number of col metadata fields (+1 for the 'id' row -- first row)
---------------------------------------------------
|id|        rhd          |          cid           |
---------------------------------------------------
|  |                     |                        |
|c |                     |                        |
|h |      (blank)        |      col_metadata      |
|d |                     |                        |
|  |                     |                        |
----------------------------------------------------
|  |                     |                        |
|r |                     |                        |
|i |    row_metadata     |          data          |
|d |                     |                        |
|  |                     |                        |
----------------------------------------------------
"""
import logging
import setup_GCToo_logger as setup_logger
import pandas as pd
import os
import GCToo

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"
logger = logging.getLogger(setup_logger.LOGGER_NAME)


def parse(file_path, nan_values=None):
    """The main method.

    Args:
        file_path: path to gct file as string
        nan_values: list of strings representing all values to be treated as NaN
            (default below)
    Returns:
        gctoo_obj: GCToo object
    """
    # Use default nan values if none given
    default_nan_values = ["#N/A", "N/A", "NA", "#NA", "NULL", "NaN", "-NaN", "nan",
                          "-nan", "#N/A!", "na", "NA", "None"]
    if nan_values is None:
        nan_values = default_nan_values

    # Read version and dimensions
    (version, num_data_rows, num_data_cols,
     num_row_metadata, num_col_metadata) = read_version_and_dims(file_path)

    # Read in metadata and data
    (row_metadata, col_metadata, data) = parse_into_3_df(file_path,
                                                         num_data_rows, num_data_cols,
                                                         num_row_metadata, num_col_metadata,
                                                         nan_values)

    # Create the gctoo object and assemble 3 component dataframes
    gctoo_obj = create_gctoo_obj(file_path, version, row_metadata, col_metadata, data)
    return gctoo_obj


def read_version_and_dims(file_path):
    # Verify that the gct path exists
    if not os.path.exists(file_path):
        err_msg = "The given path to the gct file cannot be found. gct_path: {}"
        logger.error(err_msg.format(file_path))
        raise(Exception(err_msg.format(file_path)))

    logger.info("Reading GCT: {}".format(file_path))

    # Open file
    f = open(file_path, "rb")

    # Get version from the first line
    version = f.readline().strip()

    # Check that the version is v1.3
    if version != "#1.3":
        err_msg = "Only GCT v1.3 is currently supported. version: {}"
        logger.error(err_msg.format(version))
        raise(Exception(err_msg.format(version)))

    # Read dimensions from the second line
    dims = f.readline().strip().split("\t")

    # Close file
    f.close()

    # Check that the second row is what we expect
    if len(dims) != 4:
        err_msg = "The second row of the GCT file should only have 4 entries. dims: {}"
        logger.error(err_msg.format(dims))
        raise(Exception(err_msg.format(dims)))

    # Explicitly define each dimension
    num_data_rows = int(dims[0])
    num_data_cols = int(dims[1])
    num_row_metadata = int(dims[2])
    num_col_metadata = int(dims[3])

    # Return version and dimensions
    return version, num_data_rows, num_data_cols, num_row_metadata, num_col_metadata


def parse_into_3_df(file_path, num_data_rows, num_data_cols, num_row_metadata, num_col_metadata, nan_values):
    # Read the gct file beginning with line 3
    full_df = pd.read_csv(file_path, sep="\t", skiprows=2, dtype=str,
                          na_values=nan_values, keep_default_na=False)

    # Assemble metadata dataframes
    row_metadata = assemble_row_metadata(full_df, num_col_metadata, num_data_rows, num_row_metadata)
    col_metadata = assemble_col_metadata(full_df, num_col_metadata, num_row_metadata, num_data_cols)

    # Assemble data dataframe
    data = assemble_data(full_df, num_col_metadata, num_data_rows, num_row_metadata, num_data_cols)

    # Return 3 dataframes
    return row_metadata, col_metadata, data


def assemble_row_metadata(full_df, num_col_metadata, num_data_rows, num_row_metadata):
    row_metadata_row_inds = range(num_col_metadata, num_col_metadata + num_data_rows)
    row_metadata_col_inds = range(num_row_metadata + 1)
    row_metadata = full_df.iloc[row_metadata_row_inds, row_metadata_col_inds]
    row_metadata.set_index("id", inplace=True)
    row_metadata.index.name = "rid"
    row_metadata.columns.name = "rhd"
    return row_metadata


def assemble_col_metadata(full_df, num_col_metadata, num_row_metadata, num_data_cols):
    col_metadata_row_inds = range(num_col_metadata)
    col_metadata_col_inds = [0] + range(1 + num_row_metadata, 1 + num_row_metadata + num_data_cols)
    col_metadata = full_df.iloc[col_metadata_row_inds, col_metadata_col_inds]
    col_metadata.set_index("id", inplace=True)
    col_metadata.index.name = "chd"
    col_metadata.columns.name = "cid"
    return col_metadata


def assemble_data(full_df, num_col_metadata, num_data_rows, num_row_metadata, num_data_cols):
    data_row_inds = range(num_col_metadata, num_col_metadata + num_data_rows)
    data_col_inds = [0] + range(1 + num_row_metadata, 1 + num_row_metadata + num_data_cols)
    data = full_df.iloc[data_row_inds, data_col_inds]
    data.set_index("id", inplace=True)
    data.index.name = "rid"
    data.columns.name = "cid"
    return data


def create_gctoo_obj(file_path, version, row_metadata_df, col_metadata_df, data_df):
    # Move dataframes into GCToo object; assembly occurs within GCToo.GCToo
    gctoo_obj = GCToo.GCToo(src=file_path,
                            version=version,
                            row_metadata_df=row_metadata_df,
                            col_metadata_df=col_metadata_df,
                            data_df=data_df)
    return gctoo_obj
