"""
LOG_TRANSFORM: Take the log of each value in a numpy ndarray.
Can set the base of the logarithm (b); default = 2.

"""
import utils.setup_logger as setup_logger
import logging
import argparse
import numpy

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():
	parser = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("data", help="numpy ndarray of data to be transformed") # type-check here?
	parser.add_argument("log_base", help="base of the logarithm to be applied",
		type=int, default=2)

def log_transform(data, log_base=2):

	# Data should be a numpy.ndarray
	assert isinstance(data, numpy.ndarray), ["data must be a numpy ndarray. "
		+ "type(data): %s "] % type(data)

	# How to check for NaN / weird entries?

	transformed_data = numpy.log(data) / numpy.log(log_base)
	return transformed_data

# Work to make this cmd line callable
if __name__ == "__main__":
	args = build_parser().parse_args(sys.argv[1:])
	setup_logger.setup(verbose=args.verbose)
	logger.debug("args: {}".format(args))

	log_transform(args.data, args.log_base)