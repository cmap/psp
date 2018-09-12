import unittest
import mock
import logging
import os
import broadinstitute_psp.utils.setup_logger as setup_logger
import pour

# mock setup environment variables
API_BASE_URL= "http://API_URL"
environDict = {"API_KEY": "API_KEY",
               "API_URL":API_BASE_URL,
               "PANORAMA_USER": "test_user",
               "PANORAMA_AUTH":"test_auth"}

def get_environ_item(name):
    return environDict[name]

pour.os.environ = mock.MagicMock()
pour.os.environ.__getitem__.side_effect = get_environ_item

#set up s3 mocks
s3 = mock.Mock()
s3.Object = mock.Mock()
s3.Object.return_value.get = mock.MagicMock()

s3Dict = {"Body":True}
def s3_get_item(key):
    return s3Dict[key]
s3.Object.return_value.get.__getitem__.side_effect = s3_get_item
s3.Object.return_value.get.read = mock.Mock()

pour.boto3.resource = mock.Mock(return_value=s3)


class TestPour(unittest.TestCase):

    def test_get_panorama_request_and_parse(self):
        pour.json.loads = mock.Mock(return_value={"id":"test_id"})

        bucket_name = "test_bucket"
        current_gct_key = "psp/level4/test_plate_name_LVL4.gct"

        #happy path should return request_id
        request_id = pour.get_panorama_request_and_parse(s3, bucket_name, current_gct_key)
        self.assertEqual("test_id", request_id)

    def test_get_api_entry_from_proteomics_clue(self):
        pass

    # post_update_to_proteomics_clue fully tested in utils.test_lambda_utils

if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()