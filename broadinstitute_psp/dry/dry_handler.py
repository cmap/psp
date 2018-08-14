import boto3
import botocore
import argparse
import sys
import broadinstitute_psp.utils.lambda_utils as utils
import broadinstitute_psp.dry.dry as dry
import broadinstitute_psp.utils.config_converter as config_converter

FILE_EXTENSION = ".gct"
LOCAL_LEVEL_3_GCT_NAME = "level3"
DRY_OUT_NAME = LOCAL_LEVEL_3_GCT_NAME + dry.DEFAULT_GCT_SUFFIX
LOCAL_LEVEL_2_GCT_NAME = "level2.gct"
LEVEL_3_API_SUFFIX = "/level3"
CONFIG_API_SUFFIX = "/configObj"

def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required arg
    parser.add_argument("--bucket_name", "-b", required=True,
                        help="s3 bucket where level 2 GCT is located")
    # Required arg
    parser.add_argument("--file_key", "-key", required=True,
                        help="s3 key to level 2 GCT file")
    # Required arg
    parser.add_argument("--config_dir", "-dir", required=True,
                        help="directory when psp_production config is located")

    # Required arg
    parser.add_argument("--plate_api_id", "-api", required=True,
                        help="plate id used to update API entry")

    parser.add_argument("--plate_name", "-pn", required=True,
                        help="name of plate to be processed")


    return parser

def call_dry(args):
    s3 = boto3.resource('s3')
    local_gct_path = args.config_dir + "/" + LOCAL_LEVEL_2_GCT_NAME

    download_gct_from_s3(s3, args, local_gct_path)

    try:
        config_path = check_gct_for_custom_parameters_and_set_config_path(args, local_gct_path)

    except Exception as error:
        level_3_message = "dry error: {}".format(error)
        print level_3_message
        payload = {"s3": {"message": level_3_message}}
        utils.post_update_to_proteomics_clue(LEVEL_3_API_SUFFIX, args.plate_api_id, payload)
        raise Exception(error)

    dry_args = dry.build_parser().parse_args(["-i", local_gct_path, "-p", config_path, "-o", args.config_dir, "-ob", LOCAL_LEVEL_3_GCT_NAME])
    (level_3_key, config_key) = create_s3_keys(args)

    try:
        level_3_gct = dry.main(dry_args)
        print level_3_gct

    except Exception as error:
        level_3_message = "dry error: {}".format(error)
        print level_3_message
        payload = {"s3": {"message": level_3_message}}
        utils.post_update_to_proteomics_clue(LEVEL_3_API_SUFFIX, args.plate_api_id, payload)
        raise Exception(error)

    gct_location = args.config_dir + "/" + DRY_OUT_NAME
    upload_file_to_s3(s3, args, gct_location, level_3_key, LEVEL_3_API_SUFFIX)

    # todo: if first upload fails, will not proceed to config upload, maybe this is fine?
    upload_file_to_s3(s3, args, config_path, config_key, CONFIG_API_SUFFIX)

    dry_success_payload = {"status": "created LVL 3 GCT"}
    utils.post_update_to_proteomics_clue("", args.plate_api_id, dry_success_payload)
    return "Success!"

def download_gct_from_s3(s3, args, local_level_2_path):
    try:
        print 'Reading file {} from bucket {}'.format(args.file_key, args.bucket_name)
        s3.Bucket(args.bucket_name).download_file(args.file_key, local_level_2_path)

    except botocore.exceptions.ClientError as e:

        if e.response['Error']['Code'] == "404":
            level_3_message = "The LVL2 GCT located at {} from bucket {} does not exist".format(args.file_key, args.bucket_name)
            print level_3_message

        else:
            level_3_message = "failed to download LVL2 GCT located at {} from bucket {}".format(args.file_key, args.bucket_name)
            print level_3_message

        payload = {"s3": {"message": level_3_message}}
        utils.post_update_to_proteomics_clue(LEVEL_3_API_SUFFIX, args.plate_api_id, payload)
        raise Exception(e)


def check_gct_for_custom_parameters_and_set_config_path(args, local_gct_path):
    assay = args.plate_name.split("_")[1].lower()
    config_path = args.config_dir + "/" + args.plate_name + ".cfg"

    diff_params = config_converter.convert_gct_to_config(assay, local_gct_path, config_path)
    if diff_params is None:
        config_path = args.config_dir + "/psp_production.cfg"

    return config_path


def create_s3_keys(args):
    filename = args.plate_name + "_LVL3" + FILE_EXTENSION

    top_level = args.file_key.split("/", 1)[0]
    #split keeps only top level directory
    level_3_key =  top_level + "/level3/" + filename
    config_key = top_level + "/config/" + args.plate_name + ".cfg"

    return level_3_key, config_key

def upload_file_to_s3(s3, args, local_path, s3_key, api_suffix):
    try:
        gct = open(local_path, 'rb')
        s3.Bucket(args.bucket_name).put_object(Key=s3_key, Body=gct)

    except boto3.exceptions.S3UploadFailedError as error:
        error_message = "s3 upload error"
        print error_message
        payload = {"s3": {"message": error_message}}
        utils.post_update_to_proteomics_clue(api_suffix, args.plate_api_id, payload)
        raise Exception(error)

    s3_url = "s3://" + args.bucket_name + "/" + s3_key
    success_payload = {"s3": {"url": s3_url}}
    utils.post_update_to_proteomics_clue(api_suffix, args.plate_api_id, success_payload)

if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    call_dry(args)