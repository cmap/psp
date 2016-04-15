import pandas as pd

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
"""


class GCToo(object):
    """Class representing parsed gct(x) objects as pandas dataframes.

    Contains 3 component dataframes (row_metadata_df, column_metadata_df,
    and data_df) as well as an assembly of these 3 into a multi index df
    that provides an alternate way of selecting data.
    """
    def __init__(self, src=None, version=None,
                 row_metadata_df=None, col_metadata_df=None, data_df=None):
        self.src = src
        self.version = version
        self.row_metadata_df = row_metadata_df
        self.col_metadata_df = col_metadata_df
        self.data_df = data_df
        self.multi_index_df = None

        # Can only assemble if three component dataframes exist
        if ((self.row_metadata_df is not None) and
                (self.col_metadata_df is not None) and
                (self.data_df is not None)):
            self.multi_index_df = self.assemble_multi_index_df()


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
        # Convert row_metadata to the row multi index. Need to transpose the metadata numpy array.
        row_multi_index = pd.MultiIndex.from_arrays(self.row_metadata_df.T.values, names=self.row_metadata_df.columns)

        # Convert col_metadata to the col multi index
        # N.B. Column metadata is transposed.
        transposed_col_metadata = self.col_metadata_df.T
        col_multi_index = pd.MultiIndex.from_arrays(transposed_col_metadata.values, names=transposed_col_metadata.index)

        # Create multi index dataframe using the values of data_df and the indexes created above
        multi_index_df = pd.DataFrame(data=self.data_df.values, index=row_multi_index, columns=col_multi_index)

        self.multi_index_df = multi_index_df
        return multi_index_df


def slice(gctoo, row_bool=None, col_bool=None):
    """Extract a subset of data from a GCToo object in a variety of ways.

    Args:
        gctoo (GCToo object)
        row_id (list of strings): if empty, will use all rid
        row_bool (list of bools): length must equal gctoo.data_df.shape[0]
        col_id (list of strings): if empty, will use all cid
        col_bool (list of bools): length must equal gctoo.data_df.shape[1]

        NOT YET IMPLEMENTED:
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








