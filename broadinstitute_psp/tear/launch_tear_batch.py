import boto3
import json
import os


def lambda_handler(event, context):
    s3 = boto3.resource('s3')
    bucket_name = event['Records'][0]['s3']['bucket']['name']
    file_key = event['Records'][0]['s3']['object']['key']
    (request_id, plate_name) = get_panorama_request_and_parse(s3, bucket_name, file_key)
    send_to_Batch(bucket_name, file_key, request_id, plate_name)


def get_panorama_request_and_parse(s3, bucket_name, current_gct_key):
    """ EXPECTS current_gct_key format: 'psp/levelX/FILE_NAME
                file_name format : 'PLATE_NAME_LVLX'
    """
    s3_dir = current_gct_key.split("/", 1)[0] + "/level2"
    gct_file_name = current_gct_key.rsplit("/", 1)[1]

    plate_name = gct_file_name.rsplit("_", 1)[0]
    panorama_file_key = s3_dir + "/" + plate_name + ".json"

    try:
        s3obj = s3.Object(bucket_name, panorama_file_key)
        panorama_file_content = s3obj.get()['Body'].read()
    except Exception as error:
        print error
        raise Exception(error)

    panorama_json = json.loads(panorama_file_content)
    request_id = panorama_json["id"]

    return (request_id, plate_name)

def send_to_Batch(bucket_name, file_key, request_id, plate_name):
    client = boto3.client('batch')
    job = plate_name + " tear"
    res = client.submit_job(
        jobName=job,
        jobQueue='Macchiato-misc',
        dependsOn=[],
        jobDefinition='psp-tear',
        containerOverrides={
            'command': [
                '/cmap/bin/tear', '--bucket_name', bucket_name, '--file_key', file_key,
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
