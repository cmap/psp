import unittest
import mock
import logging
import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.utils.lambda_utils as lambda_utils

logger = logging.getLogger(setup_logger.LOGGER_NAME)

API_BASE_URL = "http://API_URL"

OG_requests = lambda_utils.requests.put
OG_os_environ = lambda_utils.os.environ

class TestLambdaUtils(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        # mock setup environment variables
        environDict = {"API_KEY": "API_KEY", "API_URL": API_BASE_URL}

        def get_environ_item(name):
            return environDict[name]

        lambda_utils.os.environ = mock.MagicMock()
        lambda_utils.os.environ.__getitem__.side_effect = get_environ_item

        # mock setup requests
        lambda_utils.requests.put = mock.Mock()
        lambda_utils.requests.put.return_value.ok.return_value = True

    @classmethod
    def tearDownClass(cls):
        lambda_utils.requests.put = OG_requests
        lambda_utils.os.environ = OG_os_environ

    # @mock.patch("broadinstitute.lambda_utils.requests.put")
    def test_post_update_to_proteomics_clue(self):#, mock_request_put):
        # mock_request_put.return_value.ok.return_value = True

        test_id = "test_id"
        r = lambda_utils.post_update_to_proteomics_clue("/suffix", test_id,  {"payload":"this"})


        #todo: NoneType is not iterable
        args, kwargs = lambda_utils.requests.put.call_args
        self.assertEqual(args[0], API_BASE_URL + "/" + test_id + "/suffix")
        self.assertEqual(kwargs["json"], {"payload":"this"})
        self.assertEqual(kwargs["headers"], {"user_key":"API_KEY"})

        # Does not test response object
        # lambda_utils.requests.put = for_resetting


if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()
