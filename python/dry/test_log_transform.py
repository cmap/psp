import unittest
import logging
import utils.setup_logger as setup_logger
import numpy as np
import log_transform

logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestLogTransform(unittest.TestCase):
	def test_oneD_array(self):
		oneD_array = np.array([2,4,np.nan,16])

		expected_output = np.array([1,2,np.nan,4])
		actual_output = log_transform.log_transform(oneD_array, 2)

		expected_close_to_actual = np.all(np.isclose(expected_output, actual_output,
			equal_nan=True))
		self.assertTrue(expected_close_to_actual)

	def test_twoD_array(self):
		twoD_array = np.array([[2,4],[8,16]])

		expected_output = np.array([[1,2],[3,4]])
		actual_output = log_transform.log_transform(twoD_array, 2)

		expected_close_to_actual = np.all(np.isclose(expected_output, actual_output,
			equal_nan=True))
		self.assertTrue(expected_close_to_actual)

if __name__ == "__main__":
	setup_logger.setup(verbose=True)

	unittest.main()