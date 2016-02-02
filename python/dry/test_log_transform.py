import unittest
import logging
import utils.setup_logger as setup_logger
import numpy
import log_transform

logger = logging.getLogger(setup_logger.LOGGER_NAME)

class TestLogTransform(unittest.TestCase):
	def test_twoD_array(self):
		twoD_array = numpy.array([2,4,numpy.nan,16])

		expected_output = numpy.array([1,2,numpy.nan,4])
		actual_output = log_transform.log_transform(twoD_array, 2)

		expected_close_to_actual = numpy.all(numpy.isclose(expected_output, actual_output,
			equal_nan=True))
		self.assertTrue(expected_close_to_actual)

	def test_threeD_array(self):
		threeD_array = numpy.array([[2,4],[8,16]])

		expected_output = numpy.array([[1,2],[3,4]])
		actual_output = log_transform.log_transform(threeD_array, 2)

		expected_close_to_actual = numpy.all(numpy.isclose(expected_output, actual_output,
			equal_nan=True))
		self.assertTrue(expected_close_to_actual)

if __name__ == "__main__":
	setup_logger.setup(verbose=True)

	unittest.main()