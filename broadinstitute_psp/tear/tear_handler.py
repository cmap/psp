import sys
import argparse
import boto3
import botocore
import broadinstitute_psp.utils.lambda_utils as utils
import broadinstitute_psp.tear.tear as tear

FILE_EXTENSION = ".gct"
LOCAL_LEVEL_4_GCT_NAME = "level4.gct"
LOCAL_LEVEL_3_GCT_NAME = "level3.gct"
LEVEL_4_API_SUFFIX = "/level4"
CONFIG_API_SUFFIX = "/configObj"

def build_tear_handler_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required arg
    parser.add_argument("--bucket_name", "-b", required=True,
                        help="s3 bucket where level 3 GCT is located")
    # Required arg
    parser.add_argument("--file_key", "-key", required=True,
                        help="s3 key to level 3 GCT file")
    # Required arg
    parser.add_argument("--config_dir", "-dir", required=True,
                        help="directory when psp_production config is located")

    # Required arg
    parser.add_argument("--plate_api_id", "-api", required=True,
                        help="plate id used to update API entry")

    parser.add_argument("--plate_name", "-pn", required=True,
                        help="name of plate to be processed")
    return parser

def call_tear(args):
    s3 = boto3.resource('s3')

    local_config_path = args.config_dir + "/"+ args.plate_name + ".cfg"
    local_level_3_gct_path = args.config_dir + "/" + LOCAL_LEVEL_3_GCT_NAME
    local_level_4_gct_path = args.config_dir + "/" + LOCAL_LEVEL_4_GCT_NAME

    (level_4_key, config_key) = create_s3keys(args)

    download_file_from_s3(s3, args, args.file_key, local_level_3_gct_path, LEVEL_4_API_SUFFIX)
    download_file_from_s3(s3, args, config_key, local_config_path, CONFIG_API_SUFFIX)

    tear_args = tear.build_parser().parse_args(["-i", local_level_3_gct_path, "-psp_config_path", local_config_path, "-o",local_level_4_gct_path , "-v"])

    try:
        level_4_gct = tear.main(tear_args)
        print level_4_gct

    except Exception as error:
        level_4_message = "tear error: {}".format(error)
        print level_4_message
        payload = {"s3": {"message": level_4_message}}
        utils.post_update_to_proteomics_clue(LEVEL_4_API_SUFFIX, args.plate_api_id, payload)
        raise Exception(error)

    try:
        gct = open(local_level_4_gct_path, 'rb')
        s3.Bucket(args.bucket_name).put_object(Key=level_4_key, Body=gct)

    except boto3.exceptions.S3UploadFailedError as error:
        level_4_message = "s3 upload error"
        print level_4_message
        payload = {"s3": {"message": level_4_message}}
        utils.post_update_to_proteomics_clue(LEVEL_4_API_SUFFIX, args.plate_api_id, payload)
        raise Exception(error)

    s3_url = "s3://" + args.bucket_name + "/" + level_4_key
    success_payload = {"s3": {"url": s3_url}}
    utils.post_update_to_proteomics_clue(LEVEL_4_API_SUFFIX, args.plate_api_id, success_payload)

    tear_success_payload = {"status": "created LVL 4 GCT"}
    utils.post_update_to_proteomics_clue("", args.plate_api_id, tear_success_payload)
    return "Success!"

def download_file_from_s3(s3, args, s3key, local_path, api_suffix):
    try:
        print 'Reading file {} from bucket {}'.format(s3key, args.bucket_name)
        s3.Bucket(args.bucket_name).download_file(s3key, local_path)

    except botocore.exceptions.ClientError as e:

        if e.response['Error']['Code'] == "404":
            error_message = "Tear: The file located at {} from bucket {} does not exist".format(s3key, args.bucket_name)
            print error_message

        else:
            error_message = "Tear: failed to download file located at {} from bucket {}".format(s3key, args.bucket_name)
            print error_message

        payload = {"s3": {"message": error_message}}
        utils.post_update_to_proteomics_clue(api_suffix, args.plate_api_id, payload)
        raise Exception(e)

def create_s3keys(args):
    filename = args.plate_name + "_LVL4" + FILE_EXTENSION

    #split keeps only top level directory
    top_level = args.file_key.split("/", 1)[0]
    level_4_key = top_level + "/level4/" + filename
    config_key = top_level + "/config/" + args.plate_name + ".cfg"

    return level_4_key, config_key


if __name__ == "__main__":
    args = build_tear_handler_parser().parse_args(sys.argv[1:])
    call_tear(args)