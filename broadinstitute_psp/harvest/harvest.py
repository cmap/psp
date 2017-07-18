"""

On Panorama:
1. Admin > Go to module > Query
2. Click on targetedms folder
3. couldn't figure out how to display data...hmmph

Result of query of 'replicateannotation' with one RunId:
- dict with length 7; keys = ['rows', 'queryName', 'rowCount', 'formatVersion', 'columnModel', 'schemaName', 'metaData']
- len(dict["rows"]) = 2112, this is a list of 2112 things
- Each of those things is a dict:
    rep_query_res["rows"][0]
    {u'ReplicateId/Name': u'A01_acq_02', u'Name': u'id', u'ReplicateId/Id': 1079098, u'Value': u'GY1-33848-001A01', u'Id': 563487, u'ReplicateId/RunId/Id': 16355}

    rep_query_res["rows"][3]
    {u'ReplicateId/Name': u'A01_acq_02', u'Name': u'det_well', u'ReplicateId/Id': 1079098, u'Value': u'A1', u'Id': 563490, u'ReplicateId/RunId/Id': 16355}

- So it's like the column metadata was rolled out.

To reassemble...
- 'Name' will become the metadata header
- 'Value' will become whatever fills the metadata header
- 'ReplicateId/Name' (or 'ReplicateId/Id' will probably be better) is what
we'll use to figure out which rows (i.e. entries in the list of 2112)
belong together


Result of query of 'generalmoleculeannotation' with one RunId is similar.
- dict with length 7; keys = ['rows', 'queryName', 'rowCount', 'formatVersion', 'columnModel', 'schemaName', 'metaData']
- len(dict["rows"]) = 17113, this is a list of 17113 things
- Each of those things is a dict:
    peptide_query_res["rows"][0]
    {u'Name': u'pr_id', u'PeptideId/Id': 561546, u'PeptideId/PeptideModifiedSequence': u'T[+56.0]K[+56.0]QTAR', u'Value': u'BI10003', u'_labkeyurl_PeptideId/PeptideGroupId/Label': u'/labkey/targetedms/LINCS/GCP/showProtein.view?id=119057', u'_labkeyurl_PeptideId/Sequence': u'/labkey/targetedms/LINCS/GCP/showPeptide.view?id=561546', u'PeptideId/PeptideGroupId/RunId/File/Id': 3030, u'PeptideId/Sequence': u'TKQTAR', u'PeptideId/PeptideGroupId/Label': u'H3K4me0 - BI10003', u'Id': 438782}

    peptide_query_res["rows"][4]
    {u'Name': u'pr_gcp_base_peptide', u'PeptideId/Id': 561546, u'PeptideId/PeptideModifiedSequence': u'T[+56.0]K[+56.0]QTAR', u'Value': u'TKQTAR', u'_labkeyurl_PeptideId/PeptideGroupId/Label': u'/labkey/targetedms/LINCS/GCP/showProtein.view?id=119057', u'_labkeyurl_PeptideId/Sequence': u'/labkey/targetedms/LINCS/GCP/showPeptide.view?id=561546', u'PeptideId/PeptideGroupId/RunId/File/Id': 3030, u'PeptideId/Sequence': u'TKQTAR', u'PeptideId/PeptideGroupId/Label': u'H3K4me0 - BI10003', u'Id': 438786}

- So it's like the row metadata was rolled out.

To reassemble...
- 'Name' will become the metadata header
- 'Value' will become whatever fills the metadata header
- 'ReplicateId/Name' (or 'ReplicateId/Id' will probably be better) is what
we'll use to figure out which rows (i.e. entries in the list of 2112)
belong together

Note that in this case, the fields will be different for GCP and P100 (i.e.
different RunIds will have different metadata headers). That's okay, I guess.


We can pull the GCT from Skyline, but we can also pull annotation and data
separately and assemble it ourselves.

To pull the GCTs, we can see what Skyline files are in the 'runs' table,
parse those strings to figure out what the download URL for the unprocessed
GCT is, and download the GCT.

"""

# TODO(LL): Incomplete! Need to decide if we want to pull GCTs
# or assemble GCTs after pulling metadata.

import datetime
import logging
import labkey
import pandas as pd
import urllib
import os
import broadinstitute_psp.utils.setup_logger as setup_logger

logger = logging.getLogger(setup_logger.LOGGER_NAME)

LABKEY_SERVER = "panoramaweb.org"
PROJECT_NAME_BASE = "LINCS/"
CONTEXT_PATH = "labkey"
SCHEMA = "targetedms"
SKY_FILES_TABLE = "runs"
PREFIX_FOR_SKY_FILES = "https://panoramaweb.org/labkey/_webdav/LINCS/"
MIDDLE_STRING_FOR_SKY_FILES = "/@files/GCT/"
SUFFIX_FOR_SKY_FILES = ".sky.zip"

ROW_METADATA_TABLE = "generalmoleculeannotation"
ROW_METADATA_VIEW = "GCT_peptide_annotation"
ROW_METADATA_RUN_ID_PREFIX = "PeptideId/PeptideGroupId/RunId/File/Id"

COL_METADATA_TABLE = "replicateannotation"
COL_METADATA_VIEW = "GCT_replicate_annotation"
COL_METADATA_RUN_ID_PREFIX = "ReplicateId/RunId/Id"


def get_metadata(assay_type, run_id, run_id_prefix, table, view):
    """ Extract metadata corresponding to a particular run_id (i.e. a
    single plate).

    Args:
        assay_type (string): choices = {"GCP", "P100"}
        run_id (int)
        run_id_prefix (string): depends on the table being queried
        table (string): name of SQL table
        view (string): name of table view

    Returns:
        query_result (dict): w/ 7 keys. query_result["rows"] is a list with
            length N where each entry is one metadata value; this result
            is like the metadata df rolled out

    """
    # Create server context
    project_name = PROJECT_NAME_BASE + assay_type
    server_context = labkey.utils.create_server_context(
        LABKEY_SERVER, project_name, CONTEXT_PATH, use_ssl=True)

    # Only get results for a particular run_id
    qf = labkey.query.QueryFilter(run_id_prefix, run_id)

    # Make the SQL request
    query_result = labkey.query.select_rows(
        server_context, SCHEMA, table, filter_array=[qf],
        view_name=view)

    # Return the files that were found
    if len(query_result["rows"]) > 1:
        logger.info("result['rows'][0]: {}".format(query_result["rows"][0]))
        logger.info("Number of rows returned: " + str(query_result["rowCount"]))
    else:
        logger.error("Failed to load results from " + SCHEMA + "." + table)

    return query_result


def get_skyline_files(assay_type):
    """ For a given assay type, get back the filenames of all Skyline files of
    that assay type.

    Args:
        assay_type (string)

    Returns:
        sky_files (list of strings): each string is the name of a zipped
            .sky file

    """

    # Create server context
    project_name = PROJECT_NAME_BASE + assay_type
    server_context = labkey.utils.create_server_context(
        LABKEY_SERVER, project_name, CONTEXT_PATH, use_ssl=True)

    # Exclude Skyline files that were unsuccessfully uploaded
    qf = labkey.query.QueryFilter("Status", "")

    # Make the SQL query
    result = labkey.query.select_rows(
        server_context, SCHEMA, SKY_FILES_TABLE, filter_array=[qf])

    # Extract the Skyline files
    if len(result["rows"]) > 1:
        sky_files = [str(row["FileName"]) for row in result["rows"]]
        logger.debug("sky_files:\n{}".format(sky_files))
        logger.info("Number of rows returned: " + str(result["rowCount"]))
    else:
        logger.error("Failed to load results from " + SCHEMA + "." + SKY_FILES_TABLE)

    return sky_files


def create_sky_files_log(sky_files, dir_w_log_files):

    # Make the list into a dataframe
    sky_file_df = pd.DataFrame(sky_files)

    # Write to file with timestamp
    now_time = datetime.datetime.now()


    # TODO(LL): CONTINUE HERE!


    # actual_out_name = os.path.join(dir_w_log_files, out_name)

    pass






def create_urls_from_skyline_files(assay_type, sky_files, file_ext):

    # Create the URLS
    url_prefix = (PREFIX_FOR_SKY_FILES + assay_type + MIDDLE_STRING_FOR_SKY_FILES)
    urls = [url_prefix + sky_file.strip(SUFFIX_FOR_SKY_FILES) + file_ext for sky_file in sky_files]
    logger.debug("urls: {}".format(urls))

    return urls


def download_urls(urls, out_dir):

    for url in urls:
        # Assemble the name for where the URL should be saved
        save_name = os.path.basename(url)
        full_save_name = os.path.join(out_dir, save_name)

        # Retrieve the URL
        urllib.urlretrieve(url, full_save_name)
        logger.info("File saved to {}.".format(full_save_name))


def get_run_ids(wildcard):
    pass


def copy_unprocessed_gcts_from_panorama(assay_type, out_dir):

    skyline_files = get_skyline_files(assay_type)
    urls = create_urls_from_skyline_files(assay_type, skyline_files, ".gct")
    download_urls(urls, out_dir)


if __name__ == "__main__":
    # logger.info("main method is empty")
    setup_logger.setup(verbose=True)
    copy_unprocessed_gcts_from_panorama("GCP", "/cmap/data/proteomics/harvest/")
    copy_unprocessed_gcts_from_panorama("P100", "/cmap/data/proteomics/harvest/")
