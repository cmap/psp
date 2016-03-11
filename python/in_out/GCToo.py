import pandas as pd

__authors__ = 'Oana Enache, Lev Litichevskiy, Dave Lahr'
__email__ = 'whowantsit@broadinstitute.org'


class GCToo(object):
    """Class representing parsed gct(x) objects as pandas dataframes.

    Contains three component dataframes (row metadata, column metadata,
    and data) as well as an assembly of these three into a multi index dataframe
    that is useful for selecting data.
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

    # def __str__(self):
    # 	"""
    # 	Prints a string representation of a GCToo object.
    # 	"""
    # 	# to hold object_members
    # 	object_members_to_return = {}
    #
    # 	# get headers of row metadata and add to object_members_to_return
    # 	row_headers = list(self.row_metadata.columns.values)
    # 	object_members_to_return["row_metadata_headers"] = row_headers
    #
    # 	# get headers of column metadata and add to object_members_to_return
    # 	col_headers = list(self.col_metadata.columns.values)
    # 	object_members_to_return["col_metadata_headers"] = col_headers
    #
    # 	# add instance's object members to object_members_to_return
    # 	instance_object_members = {matrix_key: __dict__[matrix_key] for matrix_key == "data_matrix"}
    # 	object_members_to_return.update(instance_object_members)
    #
    # 	return " ".join(["{}:{}".format(k,v) for (k,v) in object_members_to_return.items()])



    ### TO-DO(lev): check for row and column duplicates
    # From online documentation: "An Index instance can only contain hashable objects".


    def assemble_multi_index_df(self):
        """
        Assembles three component dataframes into a multiindex dataframe.
        Sets the result to self.multi_index_df.

        There are pros to cons to using a multiindex df. An advantage is
        that it is easy to select data. A downside is that it is slightly
        slower than indexing using separate dfs.

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
        col_multi_index = pd.MultiIndex.from_arrays(self.col_metadata_df.values, names=self.col_metadata_df.index)

        # Create multi index dataframe using the values of data_df and the indexes created above
        multi_index_df = pd.DataFrame(data=self.data_df.values, index=row_multi_index, columns=col_multi_index)

        self.multi_index_df = multi_index_df
        return multi_index_df
