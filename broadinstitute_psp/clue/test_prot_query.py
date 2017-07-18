import logging
import os
import unittest

import broadinstitute_psp.clue.prot_query as prot_query
import broadinstitute_psp.utils.setup_logger as setup_logger

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestProtQuery(unittest.TestCase):
	def test_get_file_from_s3(self):

		out_path = "clue/ftest_read_yml_from_s3.yml"
		example_yml_on_s3 = "https://s3.amazonaws.com/data.clue.io/psp/example_user_input.yml"
		prot_query.get_file_from_s3(example_yml_on_s3, out_path)

		with open(out_path, "r") as f:
			line1 = f.readline()
			line2 = f.readline()

		self.assertEqual(line1, '[user_input]\n')
		self.assertEqual(line2, 'assay = GCP\n')

		# Clean up
		os.remove(out_path)

	def test_read_config_file(self):

		config_path = "clue/example_user_input.yml"
		[_, introspect, _, fae, out_dir, _] = prot_query.read_config_file(config_path)

		self.assertEqual(introspect, True)
		self.assertEqual(fae, ["pert_id", "cell_id", "pert_time"])
		self.assertEqual(out_dir, "527ef1c3-8859-44fb-a5b1-27f4a060a38d")

	def test_main(self):

		config_path = "https://s3.amazonaws.com/data.clue.io/psp/example_user_input.yml"
		args = prot_query.build_parser().parse_args("-u {}".format(config_path).split())
		prot_query.main(args)


if __name__ == '__main__':
	setup_logger.setup(verbose=True)
	unittest.main()
