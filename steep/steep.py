import pandas as pd
import numpy as np
import utils.psp_utils as utils

"""
Input is a merged gct file containing all data.

Output is 2 gct files: 1 for pairwise similarities and 1 for pairwise
connectivities.

"""

# Args
in_gct = "~/psp_data/all_gcp.gct"
psp_config = "~/code/broadinstitute.psp/example_psp.cfg"

# Read gct and config file
(gct, config_io, config_metadata, config_parameters) = (
    utils.read_gct_and_config_file(in_gct, psp_config))

# Compute pairwise similarity matrix


# Assemble metadata for similarity matrix

#

