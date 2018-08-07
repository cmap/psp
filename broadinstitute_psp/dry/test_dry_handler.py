import unittest
import logging
import mock
import broadinstitute_psp.utils.setup_logger as setup_logger
import dry_handler as dh
from boto3.exceptions import S3UploadFailedError
from botocore.exceptions import ClientError

logger = logging.getLogger(setup_logger.LOGGER_NAME)

dh.utils.post_update_to_proteomics_clue = mock.Mock()

# mock setup environment variables
environDict = {"API_KEY": "API_KEY", "API_URL": "API_URL"}

def get_environ_item(name):
    return environDict[name]


dh.utils.os.environ = mock.MagicMock()
dh.utils.os.environ.__getitem__.side_effect = get_environ_item

class TestDryHandler(unittest.TestCase):

    @staticmethod
    def setup_args():
        args = mock.Mock(bucket_name = "test_bucket",
                        file_key = "psp/level2/test_file_key.gct",
                        config_dir = "/this/dir",
                        plate_api_id = "test_api_id",
                        plate_name = "test_plate_name")

        return args
    def test_download_gct_from_s3(self):
        args = TestDryHandler.setup_args()
        s3 = mock.MagicMock()
        s3.Bucket("test_bucket").download_file = mock.Mock(
            side_effect=ClientError({'Error': {'Code': '404', 'Message': 'Does Not Exist'}}, 'GetObject'))

        #unhappy path ClientError ErrorCode 404
        dh.utils.post_update_to_proteomics_clue.reset_mock()
        with self.assertRaises(Exception) as context:
            dh.download_gct_from_s3(s3, args, "/this/dir/level4.gct")


        dh.utils.post_update_to_proteomics_clue.assert_called_once()
        expected_call = mock.call(dh.LEVEL_3_API_SUFFIX, args.plate_api_id,
                                  {"s3":{"message":"The LVL2 GCT located at psp/level2/test_file_key.gct from bucket test_bucket does not exist" }})
        self.assertEqual(expected_call, dh.utils.post_update_to_proteomics_clue.call_args)

        #unhappy path ClientError ErrorCode !404
        s3.Bucket("test_bucket").download_file = mock.Mock(
            side_effect=ClientError({'Error': {'Code': '500', 'Message': 'Oops'}}, 'GetObject'))
        dh.utils.post_update_to_proteomics_clue.reset_mock()

        with self.assertRaises(Exception) as context:
            dh.download_gct_from_s3(s3, args, "/this/dir/level4.gct")

        expected_call = mock.call(dh.LEVEL_3_API_SUFFIX, args.plate_api_id,
                                  {"s3": {"message": "failed to download LVL2 GCT located at psp/level2/test_file_key.gct from bucket test_bucket"}})
        self.assertEqual(expected_call, dh.utils.post_update_to_proteomics_clue.call_args)

    @mock.patch("dry_handler.config_converter")
    def test_check_gct_for_custom_parameters_and_set_config_path(self, config_converter):
        # Setup
        args = TestDryHandler.setup_args()
        local_gct_path = "/this/dir/level3.gct"

        # No differential params in GCT
        config_converter.return_value = None
        diff_params = dh.check_gct_for_custom_parameters_and_set_config_path(args, local_gct_path)
        self.assertIsNone(diff_params)



    @mock.patch("dry_handler.download_gct_from_s3")
    def test_call_dry_happy_path(self, download_gct):
        #setup
        args = TestDryHandler.setup_args()
        dh.s3 = mock.Mock()
        dh.boto3.resource = mock.Mock(return_value=dh.s3)
        dh.config_converter.convert_gct_to_config = mock.Mock(return_value="this/dir/test_plate_name.cfg")
        dh.dry.main = mock.Mock()
        dh.open = mock.Mock(return_value="opened")
        dh.s3.Bucket("test_bucket").put_object = mock.Mock()

        # happy path should call all mocks, post twice to proteomics clue
        # with s3 location and success messages
        dh.call_dry(args)
        dh.boto3.resource.assert_called_once()
        download_gct.assert_called_once()
        dh.dry.main.assert_called_once()

        # vars() turns the call_args into a dictionary, removing from Namespace()
        dry_call = vars(dh.dry.main.call_args[0][0])

        self.assertEqual(args.config_dir, dry_call["out_dir"])
        self.assertEqual("/this/dir/test_plate_name.cfg", dry_call["psp_config_path"])
        self.assertEqual("/this/dir/level2.gct",dry_call["in_gct_path"])

        open_calls = dh.open.call_args_list
        expected_open_calls = [mock.call("/this/dir/level3.dry.processed.gct", "rb"),
                               mock.call("/this/dir/test_plate_name.cfg", "rb")]
        self.assertEqual(open_calls, expected_open_calls)

        put_object_calls = dh.s3.Bucket("test_bucket").put_object.call_args_list
        expected_put_objected_calls = [mock.call(Key="psp/level3/test_plate_name_LVL3.gct", Body="opened"),
                                       mock.call(Key="psp/config/test_plate_name.cfg", Body="opened")]
        self.assertEqual(put_object_calls, expected_put_objected_calls)

        clue_posts = dh.utils.post_update_to_proteomics_clue.call_args_list
        expected_posts = [mock.call("/level3", args.plate_api_id, {"s3":{"url":"s3://test_bucket/psp/level3/test_plate_name_LVL3.gct"}}),
                          mock.call("/config", args.plate_api_id, {"s3":{"url":"s3://test_bucket/psp/config/test_plate_name.cfg"}}),
                          mock.call("", args.plate_api_id, {"status":"created LVL 3 GCT"})]
        self.assertEqual(clue_posts, expected_posts)

    @mock.patch("dry_handler.download_gct_from_s3")
    def test_call_dry_unhappy_path_dry_failure(self, download_gct):
        #setup
        args = TestDryHandler.setup_args()
        dh.s3 = mock.Mock()
        dh.boto3.resource = mock.Mock(return_value=dh.s3)
        dh.dry.main = mock.Mock(side_effect=Exception("failure"))
        dh.utils.post_update_to_proteomics_clue.reset_mock()

        # unhappy path, dry failure should call mocks, raise Exception, and post error
        with self.assertRaises(Exception) as context:
            dh.call_dry(args)
        dh.boto3.resource.assert_called_once()
        download_gct.assert_called_once()
        dh.dry.main.assert_called_once()
        self.assertEqual(str(context.exception[0]), "failure")

        dry_call = vars(dh.dry.main.call_args[0][0])

        self.assertEqual(args.config_dir, dry_call["out_dir"])
        self.assertEqual("/this/dir/test_plate_name.cfg", dry_call["psp_config_path"])
        self.assertEqual("/this/dir/level2.gct", dry_call["in_gct_path"])

        post = dh.utils.post_update_to_proteomics_clue.call_args_list[0]
        expected_post = mock.call("/level3", args.plate_api_id, {"s3":{"message":"dry error: failure"}})

        self.assertEqual(expected_post, post)

    @mock.patch("dry_handler.download_gct_from_s3")
    def test_call_dry_unhappy_path_s3_upload_failure(self, download_gct):
        # setup
        args = TestDryHandler.setup_args()
        dh.s3 = mock.Mock(name="s3 mock")
        dh.boto3.resource = mock.Mock(return_value=dh.s3)
        dh.dry.main = mock.Mock(name="dry main mock")
        dh.open = mock.Mock(name="open mock", return_value="opened")
        dh.s3.Bucket("test_bucket").put_object = mock.Mock(side_effect=S3UploadFailedError)
        dh.utils.post_update_to_proteomics_clue.reset_mock()

        with self.assertRaises(Exception) as context:
            dh.call_dry(args)
        self.assertEqual(S3UploadFailedError, type(context.exception[0]))

        dh.boto3.resource.assert_called_once()
        download_gct.assert_called_once()
        dh.dry.main.assert_called_once()

        dry_call = vars(dh.dry.main.call_args[0][0])

        self.assertEqual(args.config_dir, dry_call["out_dir"])
        self.assertEqual("/this/dir/test_plate_name.cfg", dry_call["psp_config_path"])
        self.assertEqual("/this/dir/level2.gct", dry_call["in_gct_path"])

        # exception itself is empty

        dh.open.assert_called_once()
        dh.s3.Bucket("test_bucket").put_object.assert_called_once()

        s3_put_args, s3_put_kwargs = dh.s3.Bucket("test_bucket").put_object.call_args
        print s3_put_args, s3_put_kwargs
        self.assertEqual("opened", s3_put_kwargs["Body"])
        self.assertEqual("psp/level3/test_plate_name_LVL3.gct", s3_put_kwargs["Key"])

        post = dh.utils.post_update_to_proteomics_clue.call_args_list[0]
        expected_post = mock.call("/level3", args.plate_api_id, {"s3": {"message": "s3 upload error"}})

        self.assertEqual(expected_post, post)


    def test_create_s3_keys(self):
        args = TestDryHandler.setup_args()
        returned_keys = dh.create_s3_keys(args)
        expected_keys = ("psp/level3/" + args.plate_name +"_LVL3" + dh.FILE_EXTENSION, "psp/config/"+ args.plate_name + ".cfg")
        self.assertEqual(returned_keys, expected_keys)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()