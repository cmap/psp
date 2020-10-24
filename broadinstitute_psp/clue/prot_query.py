"""
prot_query.py

This bash script consumes one configuration file (.yml) with user-defined
inputs. An example can be found under clue/example_user_input.yml. This file
lives on S3, so we need to download it.

This script is a thin wrapper for external_query_many.py. It takes care of the
user-defined config file and grabs the external, uploaded GCT from S3. A second
config file (psp_on_clue.yml) is handled by external_query_many.py.

"""
import ConfigParser
import argparse
import logging
import os
import requests
import sys
import StringIO

import broadinstitute_psp.external_query.external_query_many as eqm
import broadinstitute_psp.utils.setup_logger as setup_logger

__author__ = "Lev Litichevskiy"
__email__ = "lev@broadinstitute.org"

logger = logging.getLogger(setup_logger.LOGGER_NAME)

DUMMY_SECTION_NAME = "DUMMY_SECTION"
NAME_OF_USER_INPUT_YML = "config.yaml"


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--user_input_yml", "-u", required=True,
                        help=".yml file with user-defined inputs")
    parser.add_argument("--out_dir", "-o", default=None,
                        help=("output directory that overrides " +
                              "what's provided in user_input_yml"))
    parser.add_argument("--psp_on_clue_yml", "-p", default=None,
                        help=("path to local YML file that overrides " +
                              "what's provided in user_input_yml"))
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")
    return parser


def main(args):
    # Grab user_input_yml from S3 and store as a string
    user_input_yml_as_string =''
    if args.user_input_yml.lower().startswith("http"):
        user_input_yml_as_string = get_yml_from_s3(args.user_input_yml)
    else:
        user_input_yml_as_string = get_yml_file_local(args.user_input_yml)



    # Extract what I need from the config file string
    (assay, introspect, s3_gct_path, fae, out_dir, psp_on_clue_config_path) = (read_config_file(user_input_yml_as_string))

    # Override user-provided inputs
    if args.out_dir:
        out_dir = args.out_dir
    if args.psp_on_clue_yml:
        psp_on_clue_config_path = args.psp_on_clue_yml

    # Check if output directory exists; if not, make it
    if os.path.exists(out_dir):
        logger.info("Output directory already exists.")
    else:
        os.makedirs(out_dir)
        logger.info("Created output directory. out_dir: {}".format(out_dir))

    # Save user_input_yml into out_dir
    local_yml_path = os.path.join(out_dir, NAME_OF_USER_INPUT_YML)
    save_yml_to_file(user_input_yml_as_string, local_yml_path)

    # Download GCT from S3 and save it in out_dir
    local_gct_path = get_gct_from_s3(s3_gct_path, out_dir)

    # Call external_query_many with arguments
    eqm_args = argparse.Namespace(
        assay=assay,
        introspect=introspect,
        external_gct_path=local_gct_path,
        out_dir=out_dir,
        all=False,
        psp_on_clue_config_path=psp_on_clue_config_path,
        fields_to_aggregate_for_external_profiles=fae)

    eqm.main(eqm_args)

def get_yml_file_local(local_file):
    with open(local_file, 'r') as file:
        data = file.read()
    return data

def get_yml_from_s3(s3_path):

    # Read file in as a string
    file_as_string = requests.get(s3_path).text

    return file_as_string


def read_config_file(config_as_string):

    # Prepend a section header; ConfigParser doesn't work for files without
    # section headers
    new_top_line = "[" + DUMMY_SECTION_NAME + "]\n"
    config_as_string_with_header = new_top_line + config_as_string

    buf = StringIO.StringIO(config_as_string_with_header)
    config_parser = ConfigParser.RawConfigParser()
    config_parser.readfp(buf)

    # Unpack config file
    assay = config_parser.get(DUMMY_SECTION_NAME, "assay")
    introspect = config_parser.getboolean(DUMMY_SECTION_NAME, "introspect")
    s3_gct_path = config_parser.get(DUMMY_SECTION_NAME, "input_file")
    fae = eval(config_parser.get(DUMMY_SECTION_NAME, "fields_to_aggregate"))
    out_dir = config_parser.get(DUMMY_SECTION_NAME, "out_dir")
    psp_on_clue_config_path = config_parser.get(DUMMY_SECTION_NAME, "psp_on_clue_yml")

    return assay, introspect, s3_gct_path, fae, out_dir, psp_on_clue_config_path


def save_yml_to_file(yml_string, out_path):

    with open(out_path, "w") as f:
        f.write(yml_string)

    return None


def get_gct_from_s3(s3_path, out_dir):

    if s3_path.lower().startswith("http"):
        name_of_gct_file = s3_path.split("/")[-1]
        # Configure the output filename
        local_gct_path = os.path.join(out_dir, name_of_gct_file)

        # Read file in as a string
        file_as_string = requests.get(s3_path).text

        with open(local_gct_path, "w") as f:
            f.write(file_as_string)
    else:
        local_gct_path = s3_path
    # Get the last part of the S3 name

    return local_gct_path


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)

    main(args)
