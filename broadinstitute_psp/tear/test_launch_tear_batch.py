import unittest
import logging
import mock
import broadinstitute_psp.utils.setup_logger as setup_logger
import broadinstitute_psp.tear.launch_tear_batch as ltb
logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestLaunchDryBatch(unittest.TestCase):

    def test_get_panorama_request_and_parse(self):
        current_gct_key = "psp/level3/this_plate_name_LVL3.gct"
        bucket_name = "bucket_name"
        s3 = mock.Mock()
        s3.Object = mock.Mock()
        s3.Object.return_value.get = mock.MagicMock()

        s3Dict = {"Body":True}
        def s3_get_item(key):
            return s3Dict[key]
        s3.Object.return_value.get.__getitem__.side_effect = s3_get_item

        # happy path
        s3.Object.return_value.get.read = mock.Mock()
        ltb.json.loads = mock.Mock(return_value={"id": "fake_id"})

        returned_tuple = ltb.get_panorama_request_and_parse(s3, bucket_name, current_gct_key)
        expected_tuple = ("fake_id", "this_plate_name")
        self.assertEqual(expected_tuple, returned_tuple)

        #unhappy path s3 reading exception
        s3.Object.side_effect = Exception("failure")
        with self.assertRaises(Exception) as context:
            returned_tuple = ltb.get_panorama_request_and_parse(s3, bucket_name, current_gct_key)
        self.assertEqual("failure", context.exception[0][0])

if __name__ == "__main__":
    setup_logger.setup(verbose=True)
    unittest.main()