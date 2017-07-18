### HOWTO:

# This is a thin wrapper to call the Python script prot_query.py on Clue.

# 1. cd to the psp/broadinstitute_psp directory
# 2. Call this script.

# N.B. Need to have the wget utility function!


# This bash script runs on Clue to perform an external query. Optionally, it
# performs introspect on the external profiles. The config file
# specifies the location of required signature and similarity files.
# Before running this script, I'd advise copying clue/psp_on_clue.yml to some
# other location, changing entries as needed, and then pointing
# PSP_ON_CLUE_CONFIG_PATH to your copied config file.

#----------#

# Activate conda environment
source activate psp

python clue/prot_query.py $USER_INPUT_YML_PATH $PSP_ON_CLUE_YML_PATH

# Deactivate conda environment
source deactivate