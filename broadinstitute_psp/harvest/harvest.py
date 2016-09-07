import logging
import labkey
import utils.setup_logger as setup_logger
import urllib
import os

logger = logging.getLogger(setup_logger.LOGGER_NAME)

assay_type = "GCP"
sql_table = "runs"

labkey_server = "panoramaweb.org"
project_name_base = "LINCS/"
context_path = "labkey"
schema = "targetedms"

# TODO(lev): refactor. There should be get_metadata and get_data methods.


def get_metadata(assay_type, run_id):
    """ Extract metadata corresponding to a particular run_id (i.e. a
    single plate).

    Args:
        assay_type (string): choices = {"GCP", "P100"}
        run_id (int)

    Returns:
        query_result (list of dicts): len = #rows x #row_metadata_headers;
            keys of each dict entry correspond to metadata headers, values
            correspond to the actual metadata

    """
    # Create server context
    project_name = project_name_base + assay_type
    table = "generalmoleculeannotation"
    view = "GCT_peptide_annotation"
    server_context = labkey.utils.create_server_context(
        labkey_server, project_name, context_path, use_ssl=True)

    # Only get results for a particular run_id
    qf = labkey.query.QueryFilter("PeptideId/PeptideGroupId/RunId/File/Id",
                                  run_id)

    # Make the SQL request
    result = labkey.query.select_rows(
        server_context, schema, table, filter_array=[qf], view_name=view)

    # Return the files that were found
    if len(result["rows"]) > 1:
        logger.info("result['rows'][0]: {}".format(result["rows"][0]))
        logger.info("Number of rows returned: " + str(result["rowCount"]))
    else:
       logger.error("Failed to load results from " + schema + "." + table)


def connect_to_col_metadata(assay_type):
    """ Get list of all Skyline files on Panorama for a certain assay type.

    Args:
        assay_type (string): choices = {"GCP", "P100"}
        run_id (int)

    Returns:
        query_result (list of dicts): len = #rows x #row_metadata_headers;
            keys of each dict entry correspond to metadata headers, values
            correspond to the actual metadata

    """
    # Create server context
    project_name = project_name_base + assay_type
    table = "generalmoleculeannotation"
    view = "GCT_peptide_annotation"
    server_context = labkey.utils.create_server_context(
        labkey_server, project_name, context_path, use_ssl=True)

    # Only get results for a particular run_id
    qf = labkey.query.QueryFilter("PeptideId/PeptideGroupId/RunId/File/Id",
                                  run_id)

    # Make the SQL request
    result = labkey.query.select_rows(
        server_context, schema, table, filter_array=[qf], view_name=view)

    # Return the files that were found
    if len(result["rows"]) > 1:
        logger.info("result['rows'][0]: {}".format(result["rows"][0]))
        logger.info("Number of rows returned: " + str(result["rowCount"]))
    else:
       logger.error("Failed to load results from " + schema + "." + table)

        # ReplicateId/RunId/Id


def get_skyline_files(assay_type):
    """ For a given assay type, get back the filenames of all Skyline files of
    that assay type.

    Args:
        assay_type (string)

    Returns:
        sky_files (list of strings)

    """

    # Create server context
    project_name = project_name_base + assay_type
    table = "runs"
    server_context = labkey.utils.create_server_context(
        labkey_server, project_name, context_path, use_ssl=True)

    # Exclude Skyline files that were unsuccessfully uploaded
    qf = labkey.query.QueryFilter("Status", "")

    # Make the SQL query
    result = labkey.query.select_rows(
        server_context, schema, table, filter_array=[qf])

    # Extract the Skyline files
    if len(result["rows"]) > 1:
        sky_files = [str(row["FileName"]) for row in result["rows"]]
        logger.debug("sky_files:\n{}".format(sky_files))
        logger.info("Number of rows returned: " + str(result["rowCount"]))
    else:
       logger.error("Failed to load results from " + schema + "." + table)

    return sky_files


def create_urls_from_skyline_files(assay_type, sky_files, suffix):
    # Create the URLS
    url_prefix = ("https://panoramaweb.org/labkey/_webdav/LINCS/" +
                  assay_type + "/@files/GCT/")
    urls = [url_prefix + sky_file.strip(".sky.zip") + suffix for sky_file in sky_files]
    logger.debug("urls: {}".format(urls))

    return urls

def download_urls(urls, out_dir):

    for url in urls:
        # Assemble the name for where the URL should be saved
        save_name = os.path.basename(url)
        full_save_name = os.path.join(out_dir, save_name)

        # Retrieve the URL
        urllib.urlretrieve(url, full_save_name)




def get_run_ids(wildcard):
    pass

setup_logger.setup(verbose=False)
# connect_to_row_metadata("GCP", 3030)

assay_type = "P100"

skyline_files = get_skyline_files(assay_type)
urls = create_urls_from_skyline_files(assay_type, skyline_files, ".gct")
download_urls(urls, "/cmap/data/proteomics/harvest/wget_unprocessed")