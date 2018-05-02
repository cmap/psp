"""
This script programmatically compares the output of dry to the processed gcts
pulled down from Panorama.
"""

import glob
import os
import numpy as np
import logging
import pandas as pd

import cmapPy.pandasGEXpress.parse as parse
import broadinstitute_psp.utils.setup_logger as setup_logger

# Location of gcts processed using R code
gct_loc_r_code = "/cmap/data/proteomics/produced_by_jjaffe_code/dry/wget_processed/"
gct_r_code_suffix = "*GCP*.gct"

# Location of gcts processed using dry
gct_loc_dry = "/cmap/data/proteomics/dry/2016-10-10/"

# Suffix for files created with dry
dry_suffix = ".gct.dry.processed.gct"

##########

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)
setup_logger.setup(verbose=False)

# Get R files
all_r_files = glob.glob(gct_loc_r_code + gct_r_code_suffix)

assert len(all_r_files) > 0, "No R files were found!"
logger.debug("{} R files found".format(len(all_r_files)))

# For each R file, see if we can find its complementary file created with dry
# Only keep the R files with a complement
r_files = []
dry_files = []

for r_file in all_r_files:
    dry_file_path = os.path.join(gct_loc_dry, os.path.basename(r_file).split(".processed.gct")[0] + dry_suffix)

    if os.path.exists(dry_file_path):
        dry_files.append(dry_file_path)
        r_files.append(r_file)
    else:
        logger.debug("This dry file could not be found. dry_file_path: {}".format(
            dry_file_path))

# Do the comparison
failed_data = []
failed_col_metadata = []
failed_row_metadata = []
mean_diffs = []
for dry_file, r_file in zip(dry_files, r_files):

    dry_gct = parse.parse(dry_file)
    r_gct = parse.parse(r_file)

    logger.debug("dry_file: {}".format(dry_file))
    logger.debug("r_file: {}".format(r_file))

    # Check col metadata
    try:
        if "P100" in dry_file:
            # P100 column metadata should have one more header
            assert dry_gct.col_metadata_df.shape[1] == r_gct.col_metadata_df.shape[1] + 1, (
                "dry_gct.col_metadata_df.shape: {}, r_gct.col_metadata_df.shape: {}").format(
                dry_gct.col_metadata_df.shape, r_gct.col_metadata_df.shape[1])

            pd.util.testing.assert_frame_equal(dry_gct.col_metadata_df.iloc[:, :-1], r_gct.col_metadata_df, obj="col_metadata_df")

        else:
            pd.util.testing.assert_frame_equal(dry_gct.col_metadata_df, r_gct.col_metadata_df, obj="col_metadata_df")
    except:
        failed_col_metadata.append(dry_file)
        logger.info("failed_col_metadata: {}".format(dry_file))


    # Check row metadata
    try:
        pd.util.testing.assert_frame_equal(dry_gct.row_metadata_df, r_gct.row_metadata_df, obj="row_metadata_df")

    except:
        failed_row_metadata.append(dry_file)

    # Return some information about the differences
    atol = 0.1
    abs_diffs = abs(dry_gct.data_df.subtract(r_gct.data_df, axis=1).values.flatten())
    logger.info("number of differences > {}: {}".format(atol, (abs_diffs > atol).sum()))
    mean_diff = np.nanmedian(abs_diffs)
    logger.info("mean_abs_diff: {}".format(mean_diff))
    mean_diffs = np.append(mean_diffs, mean_diff)

    try:
        # Check that 5 digits are correct
        pd.util.testing.assert_frame_equal(dry_gct.data_df, r_gct.data_df, obj="data_df_strict")

    except:
        failed_data.append(dry_file)
        logger.info("failed_data: {}".format(dry_file))


logger.info("\ntotal number of comparisons: {}".format(len(dry_files)))
logger.info("number of failed data comparisons: {}".format(len(failed_data)))
logger.info("number of failed col_metadata comparisons: {}".format(len(failed_col_metadata)))
logger.info("number of failed row_metadata comparisons: {}".format(len(failed_row_metadata)))
logger.info("median of average differences: {}".format(np.median(mean_diffs)))

### Results:
# - All GCP plates pass
# - All P100 plates fail the data_df comparison because they have a small number
#       of samples with pretty big deviations; however, the average difference
#       of differences for P100 plates is 0.0259
# - P100, PRM, plate 26 removes 1 sample by distance filtration, even though R code doesn't


# (2016-10-05) I think I meant to write P100, PRM, Plate 29, 06H; Plate 26
# doesn't exist for either P100 or GCP.

# (2016-10-10)
# All GCP and P100 plates pass, except for P100, PRM, plate 26
