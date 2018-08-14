import sys
import argparse
import re
import os
import ConfigParser
import psp_utils
import cmapPy.pandasGEXpress.parse as parse


def build_parser():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required arg
    parser.add_argument("--assay_type", "-a", type=str, required=True, choices=["gcp","p100"])

    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input", "-i", type=str, help="JSON string to be converted")
    input_group.add_argument("--input_path", "-path", help="file path to config or GCT with pr_processing_params column"
                                                           " in row_metadata_df to be converted into a config file")

    parser.add_argument("--output_path", "-o", help="file path to write output config to", default="./custom.cfg")

    action_group = parser.add_mutually_exclusive_group(required=True)
    action_group.add_argument("--config2json", "-c2j", action="store_true",
                       help="use modified config to create JSON for insertion into GCT")
    action_group.add_argument("--json2config", "-j2c", action="store_true",
                       help="use JSON to create a custom config file for arguments differing from production config "
                            "(NB: does not create file if JSON arguments match production configuration)"
                            "requires output path")
    action_group.add_argument("--gct2config", "-g2c", action="store_true",
                       help="use ~JSON~ from GCT to create a custom config file for arguments differing from "
                            "production config (see NB above) requires output path")

    return parser


def find_psp_production_config():
    current_location = os.getcwd()
    for root, dirs, files in os.walk(current_location):
        if "psp_production.cfg" in files:
            return os.path.join(root, "psp_production.cfg")

DEFAULT_CONFIG_PATH = find_psp_production_config()


def main(args):
    if args.config2json:
        if os.path.exists(args.input_path):
            output = convert_config_to_json(args.assay_type, args.input_path)
            print output

    if args.json2config:
        convert_json_to_config(args.assay_type, args.input, args.output_path)
    if args.gct2config:
        if os.path.exists(args.input_path):
            convert_gct_to_config(args.assay_type, args.input_path, args.output_path)


def convert_config_to_json(assay, config_path):
    (_, _, custom_config_parameters) = psp_utils.read_config_file(config_path)

    differential_parameters = check_custom_parameters_against_defaults(assay, custom_config_parameters)

    if differential_parameters != None:
        output = "{"
        for i, param in enumerate(differential_parameters.keys()):
            output += ("'" + param + "':" + str(custom_config_parameters[param]))
            if i + 1 != len(differential_parameters):
                output += ","
        output += "}"
        return output
    else:
        return None


def convert_json_to_config(assay_type, json_string, output_path):
    custom_parameters = create_dict_from_pseudojson(json_string)

    differential_parameters = check_custom_parameters_against_defaults(assay_type, custom_parameters, json=True)
    write_config(differential_parameters, output_path)


def convert_gct_to_config(assay_type, gct_path, output_path):
    """
    Use custom parameters embedded within a GCT to create a config file for processing

    Args -
        gct_path -  where to read in GCT for custom parameters
        assay_type (string) - assay used to specify parameters in config
        output_path - where to write output config
    Returns -
        Nothing - writes output to output path
    """
    gct = parse.parse(gct_path)

    # All rows have same parameters embedded, choose any
    try :
        custom_params = gct.row_metadata_df.loc[:,"pr_processing_params"].any()
    except Exception as error:
        print "GCT does not contain pr_processing_params field"
        return None

    if custom_params == "{}":
        print "GCT contains pr_processing_params field, but it is empty"
        return None

    custom_params = create_dict_from_pseudojson(custom_params)
    differential_parameters = check_custom_parameters_against_defaults(assay_type, custom_params, json=True)
    if differential_parameters is not None:
        write_config(differential_parameters, output_path)

    return differential_parameters


def create_dict_from_pseudojson(pseudoJSONString):
    """
    Panorama uses a SQL DB to hold values, which coerces JSON into string here referred to as pseudo-JSON
    Use pseudo-JSON embedded in Panorama GCT to create a dictionary of parameters

    Args -
    pseudoJSONString (string) - Panorama coerced string of embedded parameters

    Returns -
    dict (dictionary) - embedded parameters

    """
    list = re.split(r'[{,}]', pseudoJSONString)
    start = None

    # Removes "base_histone_normalization_mapping" which has a list as its value
    # and will fail to parse correctly. This is utilized on Panorama during GCT creation.
    # Used to populate "pr_normalization_peptide_id" field

    for i, entry in enumerate(list):
        if "[" in entry :
            start = i
        elif "]" in entry :
            end = i + 1

    if start is not None :
        del list[start:end]

    key_value_pairs = []
    dict = {}
    for token in list:
        if (token != "") & (token != '""'):
            key_value_pairs.append(token)

    for pair in key_value_pairs:
        (key, value) = pair.split(":")

        # Strip any whitespace and remove any extra quotations
        key = key.translate(None, '"').strip()
        dict[key] = value
    return dict


def populate_background_parameters(json_parameters_dict):

    (_, _, default_config_parameters) = psp_utils.read_config_file(DEFAULT_CONFIG_PATH)
    for param in default_config_parameters:
        if param not in json_parameters_dict:
            json_parameters_dict[param] = default_config_parameters[param]
    return json_parameters_dict


def check_custom_parameters_against_defaults(assay_type, custom_config_parameters, json=False):
    """
    Converts and evaluates set of custom parameters against default config. If json, maps parameters and
    populates background (default) parameters for non-specified parameters.

    NB: requires that default config remain unchanged from that on github

    Args -
        assay (string) - config specifies assay in arguments, thus must assay type must be prepended to name
        custom_config_parameters (dictionary) - custom parameters to check against default config
        json (boolean) - whether or not to populate background parameters

    Returns -
        differential_parameters (dictionary) - parameter that differ from default, None if all custom parameters
            match defaults

    """

    if json:
        custom_config_parameters = map_R_params(assay_type, custom_config_parameters)
        custom_config_parameters = populate_background_parameters(custom_config_parameters)

    (_, _, default_config_parameters) = psp_utils.read_config_file(DEFAULT_CONFIG_PATH)

    # Checks both keys and values for deep equality
    if custom_config_parameters == default_config_parameters:
        print "Your custom config does not differ from the production config. \n" \
              "Make sure that you have not changed your psp_production.cfg \n" \
              "To verify the integrity of your production config use 'git diff ../psp_production.cfg' OSX/Linux \n" \
              "'git diff ..\\psp_production.cfg' Windows"
        return None

    else:
        params_set = set(custom_config_parameters.items()) ^ set(default_config_parameters.items())
        differential_parameters = { param : custom_config_parameters[param] for param,_ in params_set }

    return differential_parameters


def map_R_params(assay_type, parameters):
    """
    Checks for R or Python parameters and sets parameters to appropriate PSP config naming.
    Prepends assay_type for usage in processing of specific plate config is created for.
    If parameters are not in R_conversion_mapping, assumes that assay_type is included in parameters.

    Args -
    assay_type (string) - assay which parameter is used for
    parameters (dictionary) - set of parameters to convert to

    Returns -
    parameters (dictionary) - mapped parameters to appropriate PSP processing param names
        returns None if input contains invalid parameters
    """
    R_conversion_mapping = {
        "samplePctCutoff"   : "sample_frac_cutoff",
        "probePctCutoff"    : "probe_frac_cutoff",
        "probeSDCutoff"     : "probe_sd_cutoff",
        "distSDcutoff"      : "dist_sd_cutoff"
    }
    output_params = {}
    for param in parameters.keys():

        #TODO : see about adding this to psp_production.cfg

        if param == "base_histone_normalization_mapping":
            continue

        # Does not have prepended assay_type
        elif param == "offset_bounds":
            key = param

        # Convert from R parameter to Py parameter
        elif param in R_conversion_mapping.keys():
            key = assay_type + "_" + R_conversion_mapping[param]

        #Convert non-assay-specific parameter to assay-specific parameter
        elif param in R_conversion_mapping.values():
            key = assay_type + "_" + param

        # Checks that parameter contains assay_type as first token and valid PSP param as second token when split
        else:
            if (assay_type == param.split("_")[0]) & (param.split("_",1)[1] in R_conversion_mapping.values()):
                key = param
            else :
                print "The following is not a valid parameter: {}".format(param)
                return None

        output_params[key.lower()] = parameters[param]

    return output_params


def write_config(custom_config_parameters, output_path):
    config_parser = ConfigParser.RawConfigParser()
    config_parser.read(os.path.expanduser(DEFAULT_CONFIG_PATH))

    for param in custom_config_parameters:
        config_parser.set("parameters", param, custom_config_parameters[param])

    with open(output_path, "w") as config_file:
        config_parser.write(config_file)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    # setup_logger.setup(verbose=args.verbose)
    # logger.debug("args: {}".format(args))

    main(args)