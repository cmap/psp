import unittest
import logging
import mock
import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.tear.tear_handler as th
from boto3.exceptions import S3UploadFailedError
from botocore.exceptions import ClientError

logger = logging.getLogger(setup_logger.LOGGER_NAME)

OG_post_updates = th.utils.post_update_to_proteomics_clue
OG_os_environ = th.utils.os.environ

class TestTearHandler(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # post_update_to_proteomics fully tested in lambda_utils
        th.utils.post_update_to_proteomics_clue = mock.Mock()

        # mock setup environment variables
        environDict = {"API_KEY": "API_KEY", "API_URL": "API_URL"}

        def get_environ_item(name):
            return environDict[name]

        th.utils.os.environ = mock.MagicMock()
        th.utils.os.environ.__getitem__.side_effect = get_environ_item

    @classmethod
    def tearDownClass(cls):
        th.utils.post_update_to_proteomics_clue = OG_post_updates
        th.utils.os.environ = OG_os_environ

    @staticmethod
    def setup_args():
        args = mock.Mock(bucket_name = "test_bucket",
                        file_key = "psp/level3/test_file_key.gct",
                        config_dir = "/this/dir",
                        plate_api_id = "test_api_id",
                        plate_name = "test_plate_name")

        return args

    def test_download_file_from_s3(self):
        # Setup, testing GCT only, but usage is same for config
        args = TestTearHandler.setup_args()
        th.utils.post_update_to_proteomics_clue.reset_mock()
        s3 = mock.MagicMock()
        s3.Bucket("test_bucket").download_file = mock.Mock(
            side_effect=ClientError({'Error': {'Code': '404', 'Message': 'Does Not Exist'}}, 'GetObject'))

        #unhappy path ClientError ErrorCode 404
        with self.assertRaises(Exception) as context:
            th.download_file_from_s3(s3, args, args.file_key, "/this/dir/level4.gct", "/level4")
        th.utils.post_update_to_proteomics_clue.assert_called_once()

        expected_call = mock.call(th.LEVEL_4_API_SUFFIX, args.plate_api_id,
                                  {"s3":{"message":"Tear: The file located at psp/level3/test_file_key.gct from bucket test_bucket does not exist" }})
        self.assertEqual(expected_call, th.utils.post_update_to_proteomics_clue.call_args)

        #unhappy path ClientError ErrorCode !404
        th.utils.post_update_to_proteomics_clue.reset_mock()
        s3.Bucket("test_bucket").download_file = mock.Mock(
            side_effect=ClientError({'Error': {'Code': '500', 'Message': 'Oops'}}, 'GetObject'))

        with self.assertRaises(Exception) as context:
            th.download_file_from_s3(s3, args, args.file_key, "/this/dir/level4.gct", "/level4")

        expected_call = mock.call(th.LEVEL_4_API_SUFFIX, args.plate_api_id,
                                  {"s3": {"message": "Tear: failed to download file located at psp/level3/test_file_key.gct from bucket test_bucket"}})
        self.assertEqual(expected_call, th.utils.post_update_to_proteomics_clue.call_args)

    @mock.patch("broadinstitute_psp.tear.tear_handler.download_file_from_s3")
    def test_call_tear_happy_path(self, download_file):
        #setup
        args = TestTearHandler.setup_args()
        th.s3 = mock.Mock()
        th.boto3.resource = mock.Mock(return_value=th.s3)
        th.tear.main = mock.Mock()
        th.open = mock.Mock(return_value="opened")
        th.s3.Bucket("test_bucket").put_object = mock.Mock()

        # happy path should call all mocks, post twice to proteomics clue
        # with s3 location and success messages
        th.call_tear(args)

        th.boto3.resource.assert_called_once()
        download_file_calls = download_file.call_args_list
        expected_download_file_calls = [
            mock.call(th.s3, args, "psp/level3/test_file_key.gct", "/this/dir/level3.gct", "/level4"),
            mock.call(th.s3, args, "psp/config/test_plate_name.cfg", "/this/dir/test_plate_name.cfg", "/configObj")]
        self.assertEqual(download_file_calls, expected_download_file_calls)

        th.tear.main.assert_called_once()

        # vars() turns the call_args into a dictionary, removing from Namespace()
        tear_call = vars(th.tear.main.call_args[0][0])

        out_name = args.config_dir + "/" + th.LOCAL_LEVEL_4_GCT_NAME
        self.assertEqual(out_name, tear_call["out_name"])
        self.assertEqual("/this/dir/test_plate_name.cfg", tear_call["psp_config_path"])
        self.assertEqual("/this/dir/level3.gct",tear_call["in_gct_path"])

        th.open.assert_called_once()
        th.s3.Bucket("test_bucket").put_object.assert_called_once()

        s3_put_args, s3_put_kwargs = th.s3.Bucket("test_bucket").put_object.call_args
        self.assertEqual("opened", s3_put_kwargs["Body"])
        self.assertEqual("psp/level4/test_plate_name_LVL4.gct", s3_put_kwargs["Key"])

        clue_posts = th.utils.post_update_to_proteomics_clue.call_args_list
        expected_posts = [mock.call(th.LEVEL_4_API_SUFFIX, args.plate_api_id, {"s3":{"url":"s3://test_bucket/psp/level4/test_plate_name_LVL4.gct"}}),
                        mock.call("", args.plate_api_id, {"status":"created LVL 4 GCT"})]
        self.assertEqual(clue_posts, expected_posts)

    @mock.patch("broadinstitute_psp.tear.tear_handler.download_file_from_s3")
    def test_call_tear_unhappy_path_tear_failure(self, download_file):
        #setup
        args = TestTearHandler.setup_args()
        th.s3 = mock.Mock()
        th.boto3.resource = mock.Mock(return_value=th.s3)
        th.tear.main = mock.Mock(side_effect=Exception("failure"))
        th.utils.post_update_to_proteomics_clue.reset_mock()

        # unhappy path, tear failure should call mocks, raise Exception, and post error
        with self.assertRaises(Exception) as context:
            th.call_tear(args)
        self.assertEqual(str(context.exception[0]), "failure")
        th.boto3.resource.assert_called_once()

        download_file_calls = download_file.call_args_list
        expected_download_file_calls = [
            mock.call(th.s3, args, "psp/level3/test_file_key.gct", "/this/dir/level3.gct", "/level4"),
            mock.call(th.s3, args, "psp/config/test_plate_name.cfg", "/this/dir/test_plate_name.cfg", "/configObj")]
        self.assertEqual(download_file_calls, expected_download_file_calls)

        th.tear.main.assert_called_once()

        tear_call = vars(th.tear.main.call_args[0][0])

        out_name = args.config_dir + "/" + th.LOCAL_LEVEL_4_GCT_NAME
        self.assertEqual(out_name, tear_call["out_name"])
        self.assertEqual("/this/dir/test_plate_name.cfg", tear_call["psp_config_path"])
        self.assertEqual("/this/dir/level3.gct", tear_call["in_gct_path"])

        post = th.utils.post_update_to_proteomics_clue.call_args_list[0]
        expected_post = mock.call("/level4", args.plate_api_id, {"s3":{"message":"tear error: failure"}})

        self.assertEqual(expected_post, post)

    @mock.patch("broadinstitute_psp.tear.tear_handler.download_file_from_s3")
    def test_call_tear_unhappy_path_s3_upload_failure(self, download_file):
        # setup
        args = TestTearHandler.setup_args()
        th.s3 = mock.Mock(name="s3 mock")
        th.boto3.resource = mock.Mock(return_value=th.s3)
        th.tear.main = mock.Mock(name="tear main mock")
        th.open = mock.Mock(name="open mock", return_value="opened")
        th.s3.Bucket("test_bucket").put_object = mock.Mock(side_effect=S3UploadFailedError)
        th.utils.post_update_to_proteomics_clue.reset_mock()

        # unhappy path, should call mocks, raise exception and post error to proteomics clue
        with self.assertRaises(Exception) as context:
            th.call_tear(args)
        # exception itself is empty
        self.assertEqual(S3UploadFailedError, type(context.exception[0]))

        th.boto3.resource.assert_called_once()
        download_file_calls = download_file.call_args_list

        expected_download_file_calls = [mock.call(th.s3, args, "psp/level3/test_file_key.gct", "/this/dir/level3.gct", "/level4"),
                                        mock.call(th.s3, args, "psp/config/test_plate_name.cfg", "/this/dir/test_plate_name.cfg", "/configObj")]
        self.assertEqual(download_file_calls, expected_download_file_calls)

        th.tear.main.assert_called_once()

        tear_call = vars(th.tear.main.call_args[0][0])

        out_name = args.config_dir + "/" + th.LOCAL_LEVEL_4_GCT_NAME
        self.assertEqual(out_name, tear_call["out_name"])
        self.assertEqual("/this/dir/test_plate_name.cfg", tear_call["psp_config_path"])
        self.assertEqual("/this/dir/level3.gct", tear_call["in_gct_path"])


        th.open.assert_called_once()
        th.s3.Bucket("test_bucket").put_object.assert_called_once()

        s3_put_args, s3_put_kwargs = th.s3.Bucket("test_bucket").put_object.call_args
        print s3_put_args, s3_put_kwargs
        self.assertEqual("opened", s3_put_kwargs["Body"])
        self.assertEqual("psp/level4/test_plate_name_LVL4.gct", s3_put_kwargs["Key"])

        post = th.utils.post_update_to_proteomics_clue.call_args_list[0]
        expected_post = mock.call(th.LEVEL_4_API_SUFFIX, args.plate_api_id, {"s3": {"message": "s3 upload error"}})

        self.assertEqual(expected_post, post)


    def test_create_s3_keys(self):
        args = TestTearHandler.setup_args()
        s3_keys = th.create_s3keys(args)
        expected_s3_keys = ("psp/level4/" + args.plate_name +"_LVL4" + th.FILE_EXTENSION,
                            "psp/config/" + args.plate_name + ".cfg")
        self.assertEqual(s3_keys, expected_s3_keys)


if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()