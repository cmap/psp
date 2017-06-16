### HOWTO:

# IMPORTANT: First cd to the psp/broadinstitute_psp directory.

# This bash script runs on Clue to perform an external query. Optionally, it
# performs introspect on the external profiles. The config file
# specifies the location of required QCNORM and SIM files. Before running this
# script, I'd advise copying clue/psp_on_clue.cfg to some other location,
# changing entries as needed, and then pointing PSP_ON_CLUE_CONFIG_PATH to
# your copied config file.

# ASSAY {GCP, P100}
# DO_INTROSPECT {True, False}
# EXTERNAL_GCT_PATH (wherever user-uploaded data is stored)
# OUT_DIR (directory in which to save output; sub-directory is automatically made by external_query_many.py)
# PSP_ON_CLUE_CONFIG_PATH (config file containing various variables; e.g. /home/psp/psp/broadinstitute_psp/psp_config.yml)
# FIELDS_TO_AGGREGATE (which fields specify how replicates should be aggregated; e.g. pert_id cell_id)

#----------#

ASSAY=$1
DO_INTROSPECT=$2
EXTERNAL_GCT_PATH=$3
OUT_DIR=$4
PSP_ON_CLUE_CONFIG_PATH=$5
FIELDS_TO_AGGREGATE="${@:6}"

# Activate conda environment
source activate psp

# Call Python script, indicating whether or not to also do introspect
if [ "$DO_INTROSPECT" == "True" ]; then
    python external_query/external_query_many.py -a $ASSAY -i -e $EXTERNAL_GCT_PATH -o $OUT_DIR -p $PSP_ON_CLUE_CONFIG_PATH -fae $FIELDS_TO_AGGREGATE
else
    python external_query/external_query_many.py -a $ASSAY -e $EXTERNAL_GCT_PATH -o $OUT_DIR -p $PSP_ON_CLUE_CONFIG_PATH -fae $FIELDS_TO_AGGREGATE
fi


# Deactivate conda environment
source deactivate
