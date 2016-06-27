import numpy as np
import pandas as pd
import logging
import utils.setup_logger as setup_logger

__authors__ = 'Oana Enache, Lev Litichevskiy, Dave Lahr'
__email__ = 'dlahr@broadinstitute.org'

"""
DATA:
-----------------------------
|  |          cid           |
-----------------------------
|  |                        |
|r |                        |
|i |          data          |
|d |                        |
|  |                        |
-----------------------------

ROW METADATA:
--------------------------
|id|        rhd          |
--------------------------
|  |                     |
|r |                     |
|i |    row_metadata     |
|d |                     |
|  |                     |
--------------------------

COLUMN METADATA:
N.B. The df is transposed from how it looks in a gct file.
---------------------
|id|      chd       |
---------------------
|  |                |
|  |                |
|  |                |
|c |                |
|i |  col_metadata  |
|d |                |
|  |                |
|  |                |
|  |                |
---------------------

N.B. rids and cids must be unique.

"""


class GCToo(object):
    """Class representing parsed gct(x) objects as pandas dataframes.

    Contains 3 component dataframes (row_metadata_df, column_metadata_df,
    and data_df) as well as an assembly of these 3 into a multi index df
    that provides an alternate way of selecting data.
    """
    def __init__(self, data_df=None, row_metadata_df=None, col_metadata_df=None,
                 src=None, version=1.3, logger_name=setup_logger.LOGGER_NAME):
        self.logger = logging.getLogger(logger_name)

        self.src = src
        self.version = version
        self.row_metadata_df = row_metadata_df
        self.col_metadata_df = col_metadata_df
        self.data_df = data_df
        self.multi_index_df = None

        # If all three component dataframes exist, first check that they are
        # consistent, and then assemble multi_index_df
        if ((self.row_metadata_df is not None) and
                (self.col_metadata_df is not None) and
                (self.data_df is not None)):
            self.check_component_dfs()
            self.assemble_multi_index_df()


    def __str__(self):
        """Prints a string representation of a GCToo object."""
        version = "GCT v{}\n".format(self.version)
        source = "src: {}\n".format(self.src)

        if self.data_df is not None:
            data = "data_df: [{} rows x {} columns]\n".format(
            self.data_df.shape[0], self.data_df.shape[1])
        else:
            data = "data_df: None\n"

        if self.row_metadata_df is not None:
            row_meta = "row_metadata_df: [{} rows x {} columns]\n".format(
            self.row_metadata_df.shape[0], self.row_metadata_df.shape[1])
        else:
            row_meta = "row_metadata_df: None\n"

        if self.col_metadata_df is not None:
            col_meta = "col_metadata_df: [{} rows x {} columns]".format(
            self.col_metadata_df.shape[0], self.col_metadata_df.shape[1])
        else:
            col_meta = "col_metadata_df: None\n"

        full_string = (version + source + data + row_meta + col_meta)
        return full_string


    def assemble_multi_index_df(self):
        """Assembles three component dataframes into a multiindex dataframe.
        Sets the result to self.multi_index_df.

        IMPORTANT: Cross-section ("xs") is the best command for selecting
        data. Be sure to use the flag "drop_level=False" with this command,
        or else the dataframe that is returned will not have the same
        metadata as the input.

        N.B. "level" means metadata header.
        N.B. "axis=1" indicates column annotations.

        Examples:
            1) Select the probe with pr_lua_id="LUA-3404":
            lua3404_df = multi_index_df.xs("LUA-3404", level="pr_lua_id", drop_level=False)

            2) Select all DMSO samples:
            DMSO_df = multi_index_df.xs("DMSO", level="pert_iname", axis=1, drop_level=False)
        """
	    #prepare row index
        self.logger.debug("Row metadata shape: {}".format(self.row_metadata_df.shape))
        row_copy = pd.DataFrame() if self.row_metadata_df.empty else self.row_metadata_df.copy()
        row_copy["rid"] = row_copy.index
        row_index = pd.MultiIndex.from_arrays(row_copy.T.values, names=row_copy.columns)

        #prepare column index
        self.logger.debug("Col metadata shape: {}".format(self.col_metadata_df.shape))
        col_copy = pd.DataFrame() if self.col_metadata_df.empty else self.col_metadata_df.copy()
        col_copy["cid"] = col_copy.index
        transposed_col_metadata = col_copy.T
        col_index = pd.MultiIndex.from_arrays(transposed_col_metadata.values, names=transposed_col_metadata.index)

        # Create multi index dataframe using the values of data_df and the indexes created above
        self.logger.debug("Data df shape: {}".format(self.data_df.shape))
        self.multi_index_df = pd.DataFrame(data=self.data_df.values, index=row_index, columns=col_index)

    def check_component_dfs(self):
        """Checks that rids are the same between data_df and row_metadata_df,
        that cids are the same between data_df and col_metadata_df. Also,
        check that rids and cids are unique."""

        # Check rid consistency
        assert np.array_equal(self.data_df.index.values, self.row_metadata_df.index.values), (
            ("The rids in data_df do not match the rids in row_metadata_df. " +
             "self.data_df.index.values: {}, self.row_metadata_df.index.values: {}").format(
                self.data_df.index.values, self.row_metadata_df.index.values))

        # Check cid consistency
        assert np.array_equal(self.data_df.columns.values, self.col_metadata_df.index.values), (
            ("The cids in data_df do not match the cids in col_metadata_df. " +
             "self.data_df.columns.values: {}, self.col_metadata_df.index.values: {}").format(
                self.data_df.columns.values, self.col_metadata_df.index.values))

        # Check rid uniqueness
        assert len(np.unique(self.data_df.index.values)) == len(self.data_df.index.values), (
            ("The rids must be unique. self.data_df.index.values:\n{}".format(
                self.data_df.index.values)))

        # Check cid uniqueness
        assert len(np.unique(self.data_df.columns.values)) == len(self.data_df.columns.values), (
            ("The cids must be unique. self.data_df.columns.values:\n{}".format(
                self.data_df.index.values)))


def slice(gctoo, row_bool=None, col_bool=None):
    """Extract a subset of data from a GCToo object in a variety of ways.

    Args:
        gctoo (GCToo object)
        row_bool (list of bools): length must equal gctoo.data_df.shape[0]
        col_bool (list of bools): length must equal gctoo.data_df.shape[1]

        NOT YET IMPLEMENTED:
        row_id (list of strings): if empty, will use all rid
        col_id (list of strings): if empty, will use all cid
        row_meta_field (list of strings)
        row_meta_values (list of strings)
        exclude_rid (bool): if true, select row ids EXCLUDING 'rid' (default: False)
        exclude_cid (bool): if true, select col ids EXCLUDING 'cid' (default: False)
        ridx (list of ints): select rows using row indices
        cidx (list of ints): select cols using col indices
        ignore_missing (bool): if true, ignore missing ids (default: False)

    Returns:
        out_gctoo (GCToo object): gctoo after slicing

    """
    # TODO(lev): should use mutually exclusive groups and argparse here

    # Initialize output
    out_gctoo = GCToo()

    # If row_bool is None, use all rids
    if row_bool is None:
        row_bool = [True] * gctoo.data_df.shape[0]
    if col_bool is None:
        col_bool = [True] * gctoo.data_df.shape[1]

    assert len(row_bool) == gctoo.data_df.shape[0], 'len(row_bool) must equal gctoo.data_df.shape[0]'
    assert len(col_bool) == gctoo.data_df.shape[1], 'len(row_bool) must equal gctoo.data_df.shape[1]'

    # Use boolean indexing
    out_gctoo.data_df = gctoo.data_df.iloc[row_bool, col_bool]
    out_gctoo.row_metadata_df = gctoo.row_metadata_df.iloc[row_bool, :]
    out_gctoo.col_metadata_df = gctoo.col_metadata_df.iloc[col_bool, :]

    return out_gctoo

