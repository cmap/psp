# INPUTS:
# - PSP_DIR
# - EXTERNAL_GCT_PATH
# - ASSAY
# - OUT_DIR
# - PSP_ON_CLUE_CONFIG_PATH

# Move to top directory of PSP so that import statements work correctly
cd PSP_DIR

# Activate conda environment
source activate psp

# Call Python script
python clue/steep_and_sip_external_on_clue.py -a ASSAY -e EXTERNAL_GCT_PATH -o OUT_DIR -p psp_on_clue.cfg

# TODO: store stdout?
# TODO: use timestamp and user name to create subdir in out_dir?

# Deactivate conda environment
source deactivate