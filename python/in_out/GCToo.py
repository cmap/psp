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
        data = "data_df: [{} rows x {} columns]\n".format(
            self.data_df.shape[0], self.data_df.shape[1])
        row_meta = "row_metadata_df: [{} rows x {} columns]\n".format(
            self.row_metadata_df.shape[0], self.row_metadata_df.shape[1])
        col_meta = "col_metadata_df: [{} rows x {} columns]".format(
            self.col_metadata_df.shape[0], self.col_metadata_df.shape[1])

        full_string = (version + source + data + row_meta + col_meta)
        return full_string


    def assemble_multi_index_df(self):
        """
        Assembles three component dataframes into a multiindex dataframe.
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
