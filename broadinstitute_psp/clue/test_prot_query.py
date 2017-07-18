import logging
import os
import shutil
import unittest

import broadinstitute_psp.clue.prot_query as prot_query
import broadinstitute_psp.utils.setup_logger as setup_logger

# Setup logger
logger = logging.getLogger(setup_logger.LOGGER_NAME)


class TestProtQuery(unittest.TestCase):
	def test_get_yml_from_s3(self):

		example_yml_on_s3 = "https://s3.amazonaws.com/data.clue.io/psp/example_user_input.yml"
		out_string = prot_query.get_yml_from_s3(example_yml_on_s3)

		self.assertEqual(out_string[:20], 'assay: GCP\nname: P10')

	def test_get_gct_from_s3(self):

		out_dir = "clue"
		s3_path = "https://s3.amazonaws.com/data.clue.io/psp/P100_MCF7_Jnk_inhibitors_6H_24H_n12x96_unzipped.gct"

		local_path = prot_query.get_gct_from_s3(s3_path, out_dir)

		self.assertEqual(local_path, "clue/P100_MCF7_Jnk_inhibitors_6H_24H_n12x96_unzipped.gct")
		self.assertTrue(os.path.exists(local_path))

		# Clean up
		os.remove(local_path)

	def test_read_config_file(self):

		config_as_string = 'assay: GCP\nname: P100 MCF7 inhibitors\nintrospect: true\nexternal_gct_path: https://s3.amazonaws.com/data.clue.io/psp/P100_MCF7_Jnk_inhibitors_6H_24H_n12x96_unzipped.gct\nfields_to_aggregate: ["pert_id", "cell_id", "pert_time"]\nout_dir: clue/527ef1c3\npsp_on_clue_yml: clue/psp_on_clue.yml\n'
		[_, introspect, _, fae, out_dir, _] = prot_query.read_config_file(config_as_string)

		self.assertEqual(introspect, True)
		self.assertEqual(fae, ["pert_id", "cell_id", "pert_time"])
		self.assertEqual(out_dir, "clue/527ef1c3")

	# Slow (~30 sec.)
	def test_main(self):

		out_dir = "clue/527ef1c3"

		# Delete out_dir if it exists already
		if os.path.exists(out_dir):
			shutil.rmtree(out_dir)

		config_path = "https://s3.amazonaws.com/data.clue.io/psp/example_user_input.yml"
		args = prot_query.build_parser().parse_args("-u {}".format(config_path).split())
		prot_query.main(args)

		# Clean up
		shutil.rmtree(out_dir)


if __name__ == '__main__':
	setup_logger.setup(verbose=True)
	unittest.main()
