"""
prot_query.py

This bash script consumes one .yml configuration file with user-defined inputs.
An example can be found under clue/user_input.yml. This file lives on S3, so
we need to download it.

his script is a thin wrapper for external_query_many.py that takes care of the
user-defined .yml file and grabs the external, uploaded GCT from S3. A second
.yml file (psp_on_clue.yml) is handled by external_query_many.py.

"""

import ConfigParser
import argparse
import logging
import os
import requests
import sys

import broadinstitute_psp.external_query.external_query_many as eqm
import broadinstitute_psp.utils.setup_logger as setup_logger

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)

TMP_USER_INPUT_YML_PATH = "/tmp/tmp_user_input.yml"


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--user_input_yml", "-u", required=True,
                        help=".yml file with user-defined inputs")
    return parser


def main(args):

    # Grab user_input_yml from S3 and save to temporary file
    get_file_from_s3(args.user_input_yml, TMP_USER_INPUT_YML_PATH)

    # Read in the temporary file and unpack the yml file
    (assay, introspect, s3_gct_path, fae, out_dir, psp_on_clue_config_path) = (
        read_config_file(TMP_USER_INPUT_YML_PATH))

    # Download GCT from S3
    name_of_gct_file = "s3_gct_path".split("/")[-1]
    local_gct_path = os.path.join(out_dir, name_of_gct_file)
    get_file_from_s3(s3_gct_path, local_gct_path)

    # Call external_query_many with arguments
    eqm_args = argparse.Namespace(
        assay=assay,
        introspect=introspect,
        external_gct_path=local_gct_path,
        out_dir=out_dir,
        psp_on_clue_config_path=psp_on_clue_config_path,
        fields_to_aggregate_for_external_profiles=fae)
    eqm.main(eqm_args)


def get_file_from_s3(s3_path, out_path):

    # Read file in as a string
    file_as_string = requests.get(s3_path).text

    # Write to temporary file
    with open(out_path, "w") as f:
        f.write(file_as_string)

    return None


def read_config_file(config_path):

    assert os.path.exists(config_path), (
        "Config file can't be found. config_path: {}".format(config_path))

    # Read config file
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(config_path)

    # Get user-defined inputs
    assay = config_parser.get("user_input", "assay")
    introspect = config_parser.getboolean("user_input", "introspect")
    s3_gct_path = config_parser.get("user_input", "external_gct_path")
    fae = eval(config_parser.get("user_input", "fields_to_aggregate"))

    # Get inputs that we generated or provided
    out_dir = config_parser.get("other", "out_dir")
    psp_on_clue_config_path = config_parser.get("other", "psp_on_clue_yml")

    return assay, introspect, s3_gct_path, fae, out_dir, psp_on_clue_config_path


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)
