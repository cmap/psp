import boto3
import json
import os

def lambda_handler(event, context):
    s3 = boto3.resource('s3')
    bucket_name = event['Records'][0]['s3']['bucket']['name']
    gct_file_key = event['Records'][0]['s3']['object']['key']
    (request_id, plate_name) = get_panorama_request_and_parse(s3, bucket_name, gct_file_key)
    send_to_Batch(bucket_name, gct_file_key, request_id, plate_name)


def send_to_Batch(bucket_name, gct_file_key, request_id, plate_name):
    client = boto3.client('batch')
    job = plate_name + " dry"
    res = client.submit_job(
        jobName=job,
        jobQueue='Macchiato-misc',
        dependsOn=[],
        jobDefinition='psp-dry',
        containerOverrides={
            'command': [
                '/cmap/bin/dry', '--bucket_name', bucket_name, '--file_key', gct_file_key,
                '--config_dir', 'broadinstitute_psp', '--plate_api_id', request_id,
                '--plate_name', plate_name
            ],
            'environment': [
                {
                    'name': 'API_KEY',
                    'value': os.environ['API_KEY']
                },
                {
                    'name': 'API_URL',
                    'value': os.environ['API_URL']
                }
            ]
        },
        retryStrategy={
            'attempts': 1
        })
    print "jobName: {} jobId: {}".format(res['jobName'], res['jobId'])


def get_panorama_request_and_parse(s3, bucket_name, current_gct_key):
    """ EXPECTS current_gct_key format : 'psp/levelX/{FILE_NAME}
                FILE_NAME format : 'PLATE_NAME_LVLX'
    """

    s3_dir = current_gct_key.rsplit("/", 1)[0]
    gct_file_name = current_gct_key.rsplit("/", 1)[1]

    #remove _LVL2 token and file extension
    plate_name = gct_file_name.rsplit("_", 1)[0]

    panorama_file_key = s3_dir + "/" + plate_name  + ".json"

    try:
        s3obj = s3.Object(bucket_name, panorama_file_key)
        print "reading {} from bucket {}".format(panorama_file_key, bucket_name)
        panorama_file = s3obj.get()['Body'].read()

    except Exception as error:
        print error
        raise Exception(error)

    panorama_json = json.loads(panorama_file)

    request_id = panorama_json["id"]

    return (request_id, plate_name)

