import glob
import os
import numpy as np
import logging

import parse_gctoo as pg
import utils.setup_logger as setup_logger

"""

This script programmatically compares the output of dry to the processed gcts
pulled down from Panorama.

"""

# Location of gcts processed using R code
gct_loc_r_code = "/cmap/data/proteomics/produced_by_jjaffe_code/dry/wget_processed/"
gct_r_code_suffix = "*P100*.gct"

# Location of gcts processed using dry
gct_loc_dry = "/cmap/data/proteomics/dry/"

# Suffix for files created with dry
dry_suffix = ".gct.dry.processed.gct"

##########

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)
setup_logger.setup(verbose=True)

# Get R files
all_r_files = glob.glob(gct_loc_r_code + gct_r_code_suffix)
assert len(all_r_files) > 1, "No R files were found!"
logger.info("{} R files found".format(len(all_r_files)))

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
failed_comparisons = []
for dry_file, r_file in zip(dry_files, r_files):

    dry_gct = pg.parse(dry_file)
    r_gct = pg.parse(r_file)

    logger.info("dry_file: {}".format(dry_file))
    logger.info("r_file: {}".format(r_file))

    assert dry_gct.data_df.shape == r_gct.data_df.shape, "Shape of data_df is wrong!"

    try:

        # Check data in 2 stages
        strict_atol = 0.1
        lenient_atol = 1
        try:
            # More strict
            assert np.allclose(dry_gct.data_df, r_gct.data_df, atol=strict_atol, equal_nan=True), (
                "data_df does not match up.")

        except:

            # Return number of differences above strict_atol
            differences = dry_gct.data_df.subtract(r_gct.data_df).values.flatten()
            logger.warning("number of differences > {}: {}".format(strict_atol, (differences > strict_atol).sum()))

            # Less strict (absolute tolerance of 1)
            assert np.allclose(dry_gct.data_df, r_gct.data_df, atol=lenient_atol, equal_nan=True), (
                "data_df STILL does not match up!")


        # Check col metadata
        if "P100" in dry_file:
            # P100 column metadata should have one more header
            assert dry_gct.col_metadata_df.shape[1] == r_gct.col_metadata_df.shape[1] + 1, (
                "dry_gct.col_metadata_df.shape: {}, r_gct.col_metadata_df.shape: {}").format(
                dry_gct.col_metadata_df.shape, r_gct.col_metadata_df.shape[1])

            assert dry_gct.col_metadata_df.iloc[:, :-1].equals(r_gct.col_metadata_df), (
                "col_metadata_df does not match up.")

        else:
            assert dry_gct.col_metadata_df.equals(r_gct.col_metadata_df), (
                "col_metadata_df does not match up.")

        # Check row metadata
        assert dry_gct.row_metadata_df.equals(r_gct.row_metadata_df), (
            "row_metadata_df does not match up.")

    except:
        failed_comparisons.append(dry_file)

# Return summary result
if len(failed_comparisons) == 0:
    logger.info("Yay! All files matched.")

else:
    logger.warning("The following {} files did not match up.".format(len(failed_comparisons)))
    for fail in failed_comparisons:
        print fail