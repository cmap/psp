import boto3
import botocore
import urllib
import json
import requests
import os

FILE_EXTENSION = ".gct"

def handler(event, context):
    """ Function called by lambda upon put of *.json in /psp/level2 bucket
    Reads panorama request in, pulls level2 GCT from panorama using link
    location provided to API, and saves GCT to /psp/level2 bucket
    """

    bucket_name = event['Records'][0]['s3']['bucket']['name']
    file_key = event['Records'][0]['s3']['object']['key']

    panorama_request = read_panorama_request_from_s3(bucket_name, file_key)
    print "panorama request: {}".format(panorama_request)

    # Obtain request id for updating API and create level2 GCT key
    (id, gct_s3key) = extract_data_from_panorama_request(panorama_request, file_key)

    harvest(panorama_request["level 2"]["panorama"]["url"], id, bucket_name, gct_s3key, "/level2")
    return "Success!"

def read_panorama_request_from_s3(bucket_name, file_key):
    """ Reads panorama request from s3, loads JSON,
    and returns request object
    """

    s3 = boto3.resource('s3')
    try:
        print 'Reading file {} from bucket {}'.format(file_key, bucket_name)
        panorama_file = s3.Object(bucket_name, file_key)
        file_content = panorama_file.get()['Body'].read()
    except Exception as error:
        print "HARVEST : error {} reading file {} from bucket {}".format(error, file_key, bucket_name)
        raise Exception(error)

    panorama_request = json.loads(file_content)
    return panorama_request


def harvest(panorama_url, api_id, bucket, s3key, api_suffix):
    """ Reads url location of file on Panorama, opens url, and uploads file
    to s3. Called for config and Level 2 GCT

    Args -
        panorama_url (URL) - url to pull file from
        api_id (string) - api id for posting status updates
        bucket (string) - bucket location on s3 to upload GCT
        s3key (string)  - key location on s3 to upload GCT
        api_suffix (string) - location to put api update with leading '/'
    """
    s3 = boto3.client('s3')

    # Read file on Panorama
    try:
        file = urllib.urlopen(panorama_url)

    except Exception as error:
        level_2_message = "error: {}".format(error)
        print level_2_message
        payload = {"s3": {"message": level_2_message}}
        post_update_to_proteomics_clue(api_suffix, api_id, payload)
        raise Exception(error)

    # Upload file to s3
    try:
        s3.upload_fileobj(file, Bucket=bucket, Key=s3key)

    except boto3.exceptions.S3UploadFailedError as error:
        level_2_message = "s3 upload error: {}".format(error)
        print level_2_message
        payload = {"s3": {"message": level_2_message}}
        post_update_to_proteomics_clue(api_suffix, api_id, payload)
        raise Exception(error)

    # Update /psp API with s3 location of file
    s3_url = "s3://" + bucket + "/" + s3key
    success_payload = {"s3": {"url": s3_url}}
    post_update_to_proteomics_clue(api_suffix, api_id, success_payload)

    # Update /psp API with success status
    if api_suffix == "/level2":
        harvest_success_payload = {"status": "created LVL2 GCT"}
        post_update_to_proteomics_clue("", api_id, harvest_success_payload)


def extract_data_from_panorama_request(panorama_request, key):
    """ Extract API id from panorama request JSON object for updates
    and create s3 key for level2 GCT based on panorama request
    provided name and current s3 key
    """

    plate_name = panorama_request["name"]
    filename = plate_name + "_LVL2" + FILE_EXTENSION

    # Create keys for level2 GCT and config, keeping main directory location
    top_level = key.rsplit("/", 1)[0]
    gct_key = top_level + "/" + filename

    request_id = panorama_request["id"]

    return (request_id, gct_key)


def post_update_to_proteomics_clue(url_suffix, id, payload):
    """ Copy of broadinstitute_psp.utils.lambda_utils.post_update_to_proteomics_clue
    This is copied to avoid having full psp repo inside of zipped harvest function
    uploaded to Amazon Lambda.
    """

    API_key = os.environ["API_KEY"]
    API_URL = os.environ["API_URL"] + "/" + id  + url_suffix

    headers = {'user_key': API_key}

    r = requests.put(API_URL,json=payload,headers=headers)
    print r.text
    if r.ok:
        print "successfully updated API at: {}".format(API_URL)
        return True
    else:
        print "failed to update API at: {} with response: {}".format(API_URL, r.text)
        return False