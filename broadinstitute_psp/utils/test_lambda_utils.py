import unittest
import mock
import logging
import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.utils.lambda_utils as lambda_utils

logger = logging.getLogger(setup_logger.LOGGER_NAME)

# mock setup environment variables
API_BASE_URL= "http://API_URL"
environDict = {"API_KEY" : "API_KEY", "API_URL" : API_BASE_URL }

def get_environ_item(name):
    return environDict[name]

lambda_utils.os.environ = mock.MagicMock()
lambda_utils.os.environ.__getitem__.side_effect = get_environ_item

# mock setup requests
lambda_utils.requests.put = mock.Mock()

class TestLambdaUtils(unittest.TestCase):

    @mock.patch("lambda_utils.requests.put")
    def test_post_update_to_proteomics_clue(self, mock_request_put):
        mock_request_put.return_value.ok.return_value = True
        test_id = "test_id"
        r = lambda_utils.post_update_to_proteomics_clue("/suffix", test_id,  {"payload":"this"})

        args, kwargs = lambda_utils.requests.put.call_args
        self.assertEqual(args[0], API_BASE_URL + "/" + test_id + "/suffix")
        self.assertEqual(kwargs["json"], {"payload":"this"})
        self.assertEqual(kwargs["headers"], {"user_key":"API_KEY"})

        # Does not test response object

if __name__ == "__main__":
    setup_logger.setup(verbose=True)

    unittest.main()