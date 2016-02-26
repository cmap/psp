"""
LOG_TRANSFORM: Take the log of each value in a np ndarray.
Can set the base of the logarithm (b); default = 2.

Will probably want this to be a sub-fcn of bigger thing.

"""
import utils.setup_logger as setup_logger
import logging
import argparse
import numpy as np

logger = logging.getLogger(setup_logger.LOGGER_NAME)

def build_parser():
	parser = argparse.ArgumentParser(description=__doc__,
		formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("data", help="np ndarray of data to be transformed") # type-check here?
	parser.add_argument("log_base", help="base of the logarithm to be applied",
		type=int, default=2)

def log_transform(data, log_base=2):

	# Data should be a np.ndarray
	assert isinstance(data, np.ndarray), ["data must be a np ndarray. "
		+ "type(data): %s "] % type(data)

	# How to check for NaN / weird entries?
	# Also want to check if any value is zero and convert it to np.nan

	transformed_data = np.log(data) / np.log(log_base)
	return transformed_data

# Work to make this cmd line callable
if __name__ == "__main__":
	args = build_parser().parse_args(sys.argv[1:])
	setup_logger.setup(verbose=args.verbose)
	logger.debug("args: {}".format(args))

	log_transform(args.data, args.log_base)