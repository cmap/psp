### HOWTO:

# IMPORTANT: First cd to the psp/broadinstitute_psp directory.

# This bash script runs on Clue to perform an external query. The config file
# specifies the location of required QCNORM and SIM files. Before running this
# script, I'd advise copying clue/psp_on_clue.cfg to some other location,
# changing fields as necessary, and then pointing PSP_ON_CLUE_CONFIG_PATH to
# your copied config file.

# PSP_DIR (e.g. psp/broadinstitute_psp/)
# EXTERNAL_GCT_PATH
# ASSAY {GCP, P100}
# OUT_DIR
# PSP_ON_CLUE_CONFIG_PATH (e.g. /home/psp/psp/broadinstitute_psp/psp_config.yml)

#----------#

ASSAY=$1
EXTERNAL_GCT_PATH=$2
OUT_DIR=$3
PSP_ON_CLUE_CONFIG_PATH=$4

# Activate conda environment
source activate psp

# Call Python script
python steep_and_sip_external/steep_and_sip_external_on_clue.py -a ASSAY -e EXTERNAL_GCT_PATH -o OUT_DIR -p PSP_ON_CLUE_CONFIG_PATH

# TODO: store stdout?
# TODO: use timestamp and user name to create subdir in out_dir?

# Deactivate conda environment
source deactivate