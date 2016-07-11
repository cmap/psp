import os
import glob
import utils.setup_logger as setup_logger
import logging
import pandas as pd
import in_out.parse_gctoo as parse_gctoo
import in_out.GCToo as GCToo
import in_out.write_gctoo as write_gctoo


logger = logging.getLogger(setup_logger.LOGGER_NAME)

def do_it():
	files = glob.glob(os.path.expanduser("~/psp_data/LINCS*.gct"))

	# Parse gct files
        gctoos = []
        for f in files:
                gctoos.append(parse_gctoo.parse(f))

	for g in gctoos:
		# For row metadata, need to remove pr_probe_suitability_manual and pr_probe_normalization_group
		# After that, will assume that metadata for the same rid in different files is identical
		updated_row_meta_df = g.row_metadata_df.drop(["pr_probe_suitability_manual", "pr_probe_normalization_group"], axis=1)
		g.row_metadata_df = updated_row_meta_df


	# Concat the row_metadata_dfs and remove duplicate rows (probably quite slow!)
	all_row_metadata_df_dups = pd.concat([x.row_metadata_df for x in gctoos], axis=0)
	logger.info("all_row_metadata_df_dups.shape: {}".format(all_row_metadata_df_dups.shape))
	all_row_metadata_df = all_row_metadata_df_dups.drop_duplicates()
	logger.info("all_row_metadata_df.shape: {}".format(all_row_metadata_df.shape))
	logger.debug("all_row_metadata_df:\n{}".format(all_row_metadata_df))

	# Re-sort index
	all_row_metadata_df_sorted = all_row_metadata_df.sort_index()

	# Concatenate the col_metadata_dfs
	all_col_metadata_df = pd.concat([x.col_metadata_df for x in gctoos], axis=0)
	logger.info("all_col_metadata_df.shape: {}".format(all_col_metadata_df.shape))


	# Concatenate the data_dfs
	all_data_df = pd.concat([x.data_df for x in gctoos], axis=1)
	logger.info("all_data_df.shape: {}".format(all_data_df.shape))

	# Make sure df shapes are correct
	assert all_data_df.shape[0] == all_row_metadata_df.shape[0], "Number of rows is incorrect."
	assert all_data_df.shape[1] == all_col_metadata_df.shape[0], "Number of columns is incorrect."
	
        logger.info("build GCToo of all...")
        g = GCToo.GCToo(version="1.3", row_metadata_df=all_row_metadata_df_sorted, col_metadata_df=all_col_metadata_df, data_df=all_data_df)

        logger.info("write to file...")
        write_gctoo.write(g, "output.gct")

if __name__ == "__main__":
        setup_logger.setup(verbose=True)

        do_it()
