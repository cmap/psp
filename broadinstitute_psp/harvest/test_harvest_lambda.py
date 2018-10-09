import unittest
import mock
import boto3.exceptions
import logging
import broadinstitute_psp.harvest.harvest_lambda as h
import broadinstitute_psp.utils.setup_logger as setup_logger

logger = logging.getLogger(setup_logger.LOGGER_NAME)

OG_requests = h.requests.put
OG_os_environ = h.os.environ


class TestHarvestLambda(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        h.requests.put = mock.Mock()

        # mock setup environment variables
        environDict = {"API_KEY": "API_KEY", "API_URL": "API_URL"}

        def get_environ_item(name):
            return environDict[name]

        h.os.environ = mock.MagicMock()
        h.os.environ.__getitem__.side_effect = get_environ_item

    @classmethod
    def tearDownClass(cls):
        h.requests.put = OG_requests
        h.os.environ = OG_os_environ

    @staticmethod
    def setup_panorama_request():
        req_id = "5afdf5b35a5fbe51cf571650"
        name = "LINCS_P100_DIA_Plate61_annotated_minimized_2018-02-06_17-09-06"
        panorama_request = {
            "name": name,
            "assay": "P100",
            "status": "Waiting_To_Download",
            "id": req_id,
            "level2": {
                "panorama": {
                    "method": "GET",
                    "url": "https://panoramaweb.org/lincs/LINCS-DCIC/P100/runGCTReport.view?runId=32201&reportName=GCT%20File%20P100"
                }
            },
            "config": {
                "panorama": {
                    "method": "GET",
                    "url" : ""
                }
            }
        }
        return (req_id, name, panorama_request)


    def test_read_panorama_request_from_s3(self):
        args = TestHarvestLambda.setup_panorama_request()
        s3 = mock.Mock()
        h.boto3.resource = mock.Mock(return_value=s3)
        s3.Object = mock.Mock(side_effect=Exception("failure"))

        # unhappy path, should raise exception, cannot post status update (no api id)
        with self.assertRaises(Exception) as context:
            h.read_panorama_request_from_s3("fake_bucket", "fake_file_key")
        print context.exception

        self.assertEqual("failure", context.exception[0][0])

    def test_harvest_happy(self):
        (req_id, name, panorama_req) = TestHarvestLambda.setup_panorama_request()
        s3 = mock.Mock()
        s3.upload_fileobj = mock.Mock()
        h.boto3.client = mock.Mock(return_value=s3)
        h.urllib.urlopen = mock.Mock()
        h.post_update_to_proteomics_clue = mock.Mock()

        # happy path, should post twice: s3 location of lvl2, status update
        h.harvest(panorama_req["level2"]["panorama"]["url"], req_id, "fake_bucket", "psp/level2/" + name + "_LVL2.gct", "/level2")
        post_update_mock_calls = h.post_update_to_proteomics_clue.call_args_list
        expected_calls = [mock.call("/level2", req_id, {"s3":{"url":"s3://fake_bucket/psp/level2/"+ name +"_LVL2.gct" }}),
                          mock.call("", req_id, {"status":"created LVL2 GCT"})]
        self.assertEqual(expected_calls, post_update_mock_calls)

        h.post_update_to_proteomics_clue.reset_mock()
        h.harvest(panorama_req["config"]["panorama"]["url"], req_id, "fake_bucket", "psp/config/" + name + ".cfg", "/config")
        post_update_mock_calls = h.post_update_to_proteomics_clue.call_args_list
        expected_calls = [mock.call("/config", req_id, {"s3": {"url": "s3://fake_bucket/psp/config/" + name + ".cfg"}})]
        self.assertEqual(expected_calls, post_update_mock_calls)


    def test_harvest_unhappy_urlopen(self):
        #test setup
        (req_id, name, panorama_req) = TestHarvestLambda.setup_panorama_request()
        s3 = mock.Mock()
        h.boto3.client = mock.Mock(return_value=s3)
        h.urllib.urlopen = mock.Mock(side_effect=Exception("failure"))
        h.post_update_to_proteomics_clue = mock.Mock()

        #unhappy url to panorama should call to urllib.urlopen, fail, and post failure to clue
        with self.assertRaises(Exception) as context:
            h.harvest(panorama_req["level2"]["panorama"]["url"], req_id, "fake_bucket", "psp/level2/fake_panorama_key.gct", "/level2")
        self.assertEqual(str(context.exception), "failure")
        self.assertEqual(panorama_req["level2"]["panorama"]["url"], h.urllib.urlopen.call_args[0][0])

        clue_post_args = h.post_update_to_proteomics_clue.call_args[0]
        self.assertEqual("/level2", clue_post_args[0], "unhappy path, urllib Exception, post to clue does not contain API URL suffix")
        self.assertEqual(req_id, clue_post_args[1], "unhappy path, urllib Exception, post to clue does not contain request id")

        error_message = "error: " + str(context.exception)
        expected_payload = { "s3" : { "message" : error_message } }
        self.assertEqual(expected_payload, clue_post_args[2], "unhappy path, urllib Exception, post to clue does not contain message" )

    def test_harvest_unhappy_s3_upload(self):
        #setup
        (req_id, name, panorama_req) = TestHarvestLambda.setup_panorama_request()
        s3 = mock.Mock()
        s3.upload_fileobj = mock.Mock(side_effect=boto3.exceptions.S3UploadFailedError)
        h.boto3.client = mock.Mock(return_value=s3)
        h.urllib.urlopen = mock.Mock()
        h.post_update_to_proteomics_clue = mock.Mock()

        #unhappy s3 upload should post "s3 upload error"
        with self.assertRaises(Exception) as context:
            h.harvest(panorama_req["level2"]["panorama"]["url"],req_id, "fake_bucket", "psp/level2/fake_panorama_key.gct", "/level2" )
        self.assertEqual(boto3.exceptions.S3UploadFailedError, type(context.exception[0]))
        h.post_update_to_proteomics_clue.assert_called_once()

        post_update_mock_call = h.post_update_to_proteomics_clue.call_args
        expected_call = mock.call("/level2", req_id, {"s3":{"message":"s3 upload error: "}})

        self.assertEqual(expected_call, post_update_mock_call)

    def test_extract_data_from_panorama_request(self):
        #setup
        (req_id, name, panorama_req) = TestHarvestLambda.setup_panorama_request()
        key = "psp/level2/request.json"

        gct_key = "psp/level2/" + name + "_LVL2" + h.FILE_EXTENSION

        (request_id, aws_gct_key) = h.extract_data_from_panorama_request(panorama_request=panorama_req, key=key)

        self.assertEqual(req_id, request_id)
        self.assertEqual(aws_gct_key, gct_key)

    ## post_update_to_proteomics_clue not explicitly tested here, see utils.lambda_utils

if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()