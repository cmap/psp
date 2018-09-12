#!/bin/bash
# Use > 1 to consume two arguments per pass in the loop (e.g. each
# argument has a corresponding value to go with it).
# Use > 0 to consume one or more arguments per pass in the loop (e.g.
# some arguments don't have a corresponding value to go with it such
# as in the --default example).
# note: if this is set to > 0 the /etc/hosts part is not recognized ( may be a bug )

while [[ $# > 1 ]]
do
key="$1"

case $key in
    --bucket_name)
    BUCKET_NAME="$2"
    shift # past argument
    ;;
    --file_key)
    FILE_KEY="$2"
    shift # past argument
    ;;
    --config_dir)
    CONFIG_DIR="$2"
    shift # past argument
    ;;
    --plate_api_id)
    PLATE_API_ID="$2"
    shift # past argument
    ;;
    --plate_name)
    PLATE_NAME="$2"
    shift # past argument
    ;;
    --default)
    DEFAULT=YES
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done

echo BUCKET_NAME  = "${BUCKET_NAME}"
echo FILE_KEY     = "${FILE_KEY}"
echo CONFIG_DIR   = "${CONFIG_DIR}"
echo PLATE_API_ID = "${PLATE_API_ID}"
echo PLATE_NAME   = "${PLATE_NAME}"

# Activate conda environment
source activate psp

cd /cmap/psp

python setup_psp.py develop

python /cmap/psp/broadinstitute_psp/dry/dry_handler.py -b ${BUCKET_NAME} -key ${FILE_KEY} -dir ${CONFIG_DIR} -api ${PLATE_API_ID} -pn ${PLATE_NAME}

# Deactivate conda environment
source deactivate


