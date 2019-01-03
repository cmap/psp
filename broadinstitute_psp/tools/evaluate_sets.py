"""
evaluate_sets.py

Evaluate sets within a similarity matrix. Sets can be replicates corresponding
to the same perturbation, or compounds corresponding to the same mechanism of
action, or something else. Sets can be provided in a GMT file, or they can be
constructed from the metadata using the match_fields argument. If constructed
from the metadata, each sample in the dataset can correspond to only one set,
but this restriction does not apply if using set definitions from a GMT file.
Finally, metadata can be embedded into a GCT(X) file or provided separately.

For each set, the aggregated similarity, a p-value, and a q-value are returned.
The p-value is computed by comparing the aggregated similarity to a size-matched
empirical null; in particular, the p-value is the fraction of values in the null
distribution greater than the aggregated similarity. The empirical null is
constructed by picking samples at random from the similarity matrix. The q-value
is computed using the Benjamini-Hochberg correction.

--------------------------------------------------------------------------------

For the following examples, running evaluate_sets from the command line will
create the directory specified with the '-o' option, and the following files will
be written to the output directory:
    - evaluate_sets results tsv containing info like p-value, q-value, set size,
        and any subsetting fields
    - nulls gctx (multiple nulls gctx if subsetting)
    - YML settings file

For within-python examples, evaluate_sets output is in the form of a
list of SimMatResultRecord objects. A SimMatResultRecord contains a results
dataframe, nulls gctoo, and subset field combination for a given run of
evaluate_sets.

--------------------------------------------------------------------------------

EXAMPLE 1

If null sets should be selected from entire similarity matrix and sets to
evaluate should be defined from a GMT file:

COMMAND LINE:
    python evaluate_sets.py
        -i /path/to/similarity/matrix
        -m /path/to/metadata
        -o /path/to/output/dir
        -s /path/to/gmt
        -mf    # match_field with no value so provided
                 GMT sets will match on profile ids

WITHIN PYTHON (check argparse defaults for desired args values):
    # run evaluate_sets and get result record list
    sim_mat_result_record_list = evaluate_sets_on_all_sim_mats(
        ds = /path/to/similarity/matrix,             # -i from command line example
        external_metadata_path = /path/to/metadata,  # -m from command line example
        set_definitions = /path/to/gmt,              # -s from command line example
        match_fields = None,                         # -mf from command line example
        subset_fields = args.subset_fields,
        external_metadata_id_field = args.external_metadata_id_field,
        aggregation_method = args.aggregation_method,
        null_size = args.null_size,
        sets_per_chunk = args.sets_per_chunk,
        low_memory_mode = args.low_memory_mode,
        max_set_size = args.max_set_size
    )
    # because we want to evaluate from entire similarity matrix, no subsetting
    # occurred, so we will have only one result record in list.
    record = sim_mat_result_record_list[0]

    # evaluate_sets results dataframe
    record.results_df

    # nulls gctoo
    record.nulls_gct

--------------------------------------------------------------------------------

EXAMPLE 2

If null sets should be selected from entire similarity matrix and sets to evaluate
should be selected based on metadata field combinations...
E.g. A set should be all profiles of the same compound at the same dose:

COMMAND LINE:
    python evaluate_sets.py
        -i /path/to/similarity/matrix
        -m /path/to/external/metadata
        -o /path/to/output/dir
        -mf pert_id pert_idose    # sets should comprise profiles of
                                    same compound, same dose

WITHIN PYTHON (check argparse defaults for desired args values):
    # run evaluate_sets and get result record list
    sim_mat_result_record_list = evaluate_sets_on_all_sim_mats(
        ds = /path/to/similarity/matrix,             # -i from command line example
        external_metadata_path = /path/to/metadata,  # -m from command line example
        match_fields = ["pert_id", "pert_idose"],    # -mf from command line example
        set_definitions = args.set_definitions,
        subset_fields = args.subset_fields,
        external_metadata_id_field = args.external_metadata_id_field,
        aggregation_method = args.aggregation_method,
        null_size = args.null_size,
        sets_per_chunk = args.sets_per_chunk,
        low_memory_mode = args.low_memory_mode,
        max_set_size = args.max_set_size
    )
    # because we want to evaluate from entire similarity matrix, no subsetting
    # occurred, so we will have only one result record in list.
    record = sim_mat_result_record_list[0]

    # evaluate_sets results dataframe
    record.results_df

    # nulls gctoo
    record.nulls_gct

--------------------------------------------------------------------------------

EXAMPLE 3

If input similarity matrix should be subsetted (e.g. based on cell line
and timepoint) before evaluating sets so that null sets and sets to evaluate
are selected for each subsetted similarity matrix:

COMMAND LINE:
    python evaluate_sets.py
        -i /path/to/similarity/matrix
        -m /path/to/external/metadata
        -o /path/to/output/dir
        -mf pert_id pert_idose    # sets should comprise profiles of
                                    same compound, same dose
        -sf cell_id pert_itime    # evaluate sets on each combination
                                    of cell line and timepoint.
                                    For example, if your input similarity matrix
                                    had A549 and MCF7 at 6 and 24H timepoints,
                                    each combination of these subset fields
                                        - A549 at 6H
                                        - A549 at 24H
                                        - MCF7 at 6H
                                        - MCF7 at 24H
                                    would have its own evaluate_sets calculation,
                                    with combination-specific sets and null.

WITHIN PYTHON (check argparse defaults for desired args values):
    # run evaluate_sets and get result record list
    sim_mat_result_record_list = evaluate_sets_on_all_sim_mats(
        ds = /path/to/similarity/matrix,             # -i from command line example
        external_metadata_path = /path/to/metadata,  # -m from command line example
        match_fields = ["pert_id", "pert_idose"],    # -mf from command line example
        subset_fields = ["cell_id", "pert_itime"],   # -sf from command line example
        set_definitions = args.set_definitions,
        external_metadata_id_field = args.external_metadata_id_field,
        aggregation_method = args.aggregation_method,
        null_size = args.null_size,
        sets_per_chunk = args.sets_per_chunk,
        low_memory_mode = args.low_memory_mode,
        max_set_size = args.max_set_size
    )

    # Because subsetting occurred on cell lines A549 and MCF7 and timepoints
    # 6H and 24H, we expect sim_mat_result_record_list to be a list of result
    # records for the four possible combinations of cell/timepoint.

    for record in sim_mat_result_record_list:
        # this instance variable will tell you which subset combination of
        # current result record. e.g. ("A549", "6 h")
        record.subset_field_combo

        # evaluate_sets results dataframe for current subset combination
        record.results_df

        # nulls gctoo for current subset combination
        record.nulls_gct
    
    # if you would like all results collated into a single dataframe, with
    # subset combinations annotated as additional columns (e.g. cell_id, pert_itime):
    collated_results_df = collate_results_dfs(sim_mat_result_record_list,
                                              args.subset_fields)
"""

import logging
import argparse
import os
import sys
import pandas as pd
import numpy as np
import random
import statsmodels.sandbox.stats.multicomp as multicomp
import itertools

import cmapPy.pandasGEXpress.setup_GCToo_logger as setup_logger
import cmapPy.pandasGEXpress.GCToo as GCToo
import cmapPy.pandasGEXpress.parse as cpp
import cmapPy.pandasGEXpress.write_gctx as wgx
import cmapPy.set_io.gmt as gmt

__author__ = "Lev Litichevskiy, Andrew Yang"
__email__ = "lev@broadinstitute.org"

MATCH_FIELD_NAME = "match_field"
SEPARATOR = ":"
AGG_SIM_COLUMN_NAME = "agg_sim"
SET_SIZE_COLUMN_NAME = "set_size"
PVAL_COLUMN_NAME = "pval"
QVAL_COLUMN_NAME = "qval"
NAN_STR = "NaN"
TEXT_FILE_OUT_NAME = "eval_sets.tsv"
NULL_DIST_OUT_NAME = "null_dists_n{n_cols}x{n_rows}{subset_fname_suffix}.gctx"
YAML_FILE_OUT_NAME = "eval_sets.yml"

logger = logging.getLogger(setup_logger.LOGGER_NAME)


class SimMatResultRecord():
    """
    Similarity matrix result record class to group settings for a run of
    eval_sets on a particular similarity matrix with the actual results
    of that run.

    For a given evaluate_sets reproducibility result and null dist., we must
    keep track of whether it was from a subsetted similarity matrix, and
    if so, what those subsetting parameters were. We do so with this class
    in order to avoid maintaining multiple lists of objects associated by
    a common index.

    Instance Variables:
        sample_ids (list of strings): A list of ids that will be set as the rid
            and cid parameters in parse(). When subsetting the similarity
            matrix, sample_ids will be the ids to subset. When not subsetting,
            sample_ids will simply be meta_df.index to parse the whole matrix.
        subset_field_combo (tuple of strings): A combination of parameters to
            subset on.
            e.g. ("A549", "48 h")
            If not subsetting, subset_field_combo will be a tuple of an empty
            string, ("")
        subset_fname_suffix (string): A filename suffix for a nulls gct,
            derived from the subset field combo.
            e.g. "_A549_48h" for a subsetted eval_sets run
            e.g. "" for a nonsubsetted eval_sets run
        results_df (pandas df): evaluate_sets results dataframe;
            see evaluate_sets_on_all_sim_mats() documentation for detailed
            explanation
        nulls_gct (gctoo): evaluate_sets null distribution gctoo;
            see null_dict_to_gctoo() documentation for detailed explanation
    """
    def __init__(self, sample_ids, subset_field_combo):
        self.sample_ids = sample_ids
        self.subset_field_combo = subset_field_combo
        if subset_field_combo == "":
            self.subset_fname_suffix = ""
        else:
            self.subset_fname_suffix = "_" + "_".join([str(sf) for sf in subset_field_combo]).replace(" ", "")

        # results_df and nulls_gct will be None at initialization.
        # they will be set after evaluating sets.
        self.results_df = None
        self.nulls_gct = None


def build_parser():
    """Build argument parser."""
    class RawDescriptionArgumentDefaultsHelpFormatter(argparse.ArgumentDefaultsHelpFormatter,
                                                      argparse.RawDescriptionHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=RawDescriptionArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--ds", "-i", required=True,
                        help="path to similarity matrix as GCT(X) file")
    # Optional args
    parser.add_argument("--match_fields", "-mf", default=["pert_id"], const=None, nargs="*",
        help=("name of metadata field(s) to use for identifying sets; "
              "the result of this grouping should match the entries in the "
              "GMT file (if provided); if flag used but no args given, "
              "value will be set to None and ids will be used"))
    parser.add_argument("--out_dir", "-o", default="eval_sets_output",
                        help="output directory (must not already exist)")
    parser.add_argument("--external_metadata_path", "-m", default=None,
        help=("path to external metadata tsv file; if None, column metadata "
              "must be embedded in the GCT(X)"))
    parser.add_argument("--external_metadata_id_field", "-emi", default="distil_id",
        help=("name of metadata field in external metadata file whose entries "
              "match cids in the GCT(X) file"))
    parser.add_argument("--set_definitions", "-s", default=None,
        help=("path to GMT file containing set definitions; if None, "
              "sets will be generated using match_fields"))
    parser.add_argument("--aggregation_method", "-am", default="median",
                        choices=["q75", "median", "median_of_medians"],
                        help="how to aggregate the self similarities into a single number")
    parser.add_argument("--subset_fields", "-sf", default=[], nargs="+",
        help=("name of metadata field(s) to use to subset similarity matrix; "
              "evaluate sets will be run on each combination of subsetted "
              "matrices so each subsetted result is computed with a distinct "
              "subsetted null. Results are then collated into a single result "
              "with additional columns annotating these subsetted fields"))
    parser.add_argument("--verbose", "-v", action="store_true", default=False,
                        help="whether to increase the # of messages reported")

    ### Hidden arguments

    # Whether to run in low-memory mode; if so, only subsets of the GCTX will
    # be read in at a time; pointless if using a GCT file"
    parser.add_argument("--low_memory_mode", "-lm", default=False,
                        action="store_true", help=argparse.SUPPRESS)

    # For low-memory mode, how many sets to read from the GCTX at once
    parser.add_argument("--sets_per_chunk", "-spc", default=100, type=int,
                        help=argparse.SUPPRESS)

    # Number of times to randomly sample for the null distribution
    parser.add_argument("--null_size", "-ns", default=10000, type=int,
                        help=argparse.SUPPRESS)

    # For efficiency, exclude sets above a certain size (e.g. positive and
    # negative controls
    parser.add_argument("--max_set_size", "-mss", default=200, type=int,
                        help=argparse.SUPPRESS)

    return parser


def main():
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    eval_sets_main(args)


def eval_sets_main(args):
    """ Separated from main(), which takes no args, in order to make eval_sets
    a command line tool.

    """
    validate_inputs(args)

    # Compute aggregated similarities, p-values, q-values, and set sizes
    sim_mat_result_record_list = evaluate_sets_on_all_sim_mats(
        ds=args.ds,
        external_metadata_path=args.external_metadata_path,
        set_definitions=args.set_definitions,
        external_metadata_id_field=args.external_metadata_id_field,
        match_fields=args.match_fields,
        aggregation_method=args.aggregation_method,
        null_size=args.null_size,
        sets_per_chunk=args.sets_per_chunk,
        low_memory_mode=args.low_memory_mode,
        max_set_size=args.max_set_size,
        subset_fields=args.subset_fields
    )

    results_df = collate_results_dfs(sim_mat_result_record_list,
                                     args.subset_fields)

    create_output_dir(args.out_dir)

    # Save evaluate_sets results
    full_out_text_name = os.path.join(args.out_dir, TEXT_FILE_OUT_NAME)
    results_df.to_csv(full_out_text_name, sep="\t", na_rep=NAN_STR,
                      float_format="%.4f")

    # Save null distributions
    nulls_fname_list = write_null_gcts(sim_mat_result_record_list, args.out_dir)

    # Save settings evaluate_sets was run with to a yaml file
    write_yaml(args, full_out_text_name, nulls_fname_list)


def validate_inputs(args):
    """
    Ensure that desired file inputs all exist, and that the output
      directory does not exist. Also, check that the parent
      directory for the output_directory does exist.
    """
    parent_dir = os.path.normpath(os.path.join(args.out_dir, os.pardir))
    if not os.path.isdir(parent_dir):
        raise Exception(("Parent directory for desired output directory"
                        "does not exist. parent_dir: {}").format(parent_dir))

    if os.path.isdir(args.out_dir):
        raise Exception(("Output directory already exists. "
                        "out_dir: {}").format(args.out_dir))

    if not os.path.isfile(args.ds):
        raise Exception(("input similarity matrix does not exist. "
                         "args.ds: {}").format(args.ds))

    if ((args.external_metadata_path is not None) and
        (not os.path.isfile(args.external_metadata_path))):
        raise Exception(("input external metadata path does not exist. "
                         "args.external_metadata_path: {}").format(
                             args.external_metadata_path))

    if ((args.set_definitions is not None) and
        (not os.path.isfile(args.set_definitions))):
        raise Exception(("input set definitions GMT path does not exist. "
                         "args.set_definitions: {}").format(
                             args.set_definitions))


def create_output_dir(out_dir):
    """
    Create output directory just before writing result files,
      ensuring we will not overwrite anything.

    Args:
        out_dir (string): output directory, from args.out_dir
    """
    try:
        os.mkdir(out_dir)
    except OSError as e:
        msg = ("Output directory already exists, cannot override."
               "out_dir: {}").format(out_dir)
        raise Exception(msg)


def evaluate_sets_on_all_sim_mats(ds, set_definitions,
                                  external_metadata_id_field,
                                  external_metadata_path, match_fields,
                                  subset_fields, aggregation_method,
                                  null_size, sets_per_chunk,
                                  low_memory_mode, max_set_size):
    """ Top-level function for evaluating sets, with or without subsetting input GCT(X)
    Args:
        ds (string): path to input GCT(X)
        set_definitions (string): path to GMT file defining sets
        external_metadata_id_field (string): metadata field in external GCT
            file that uniquely identifies samples
        external_metadata_path (string): path to external metadata file; if None,
            metadata must come from the GCT(X)
        match_fields (list of strings): metadata field(s) to use for matching samples to sets
        subset_fields (list of strings): metadata field(s) to use for subsetting ds
        aggregation_method (string): how to aggregate the self sims
        null_size (integer): number of times to randomly sample for the null dist
        sets_per_chunk (integer): number of sets to evaluate at a time in low-
            memory mode
        low_memory_mode (bool)
        max_set_size (integer): for efficiency, exclude sets above this size

    Returns:
        sim_mat_result_record_list (list of SimMatResultRecord): list of
            similarity matrix result records, with all attributes set.
            If no subsetting, this will be a list of a single SimMatResultRecord.
    """
    # input similarity matrix must be square.
    n_cid = len(cpp.parse(ds, col_meta_only=True))
    n_rid = len(cpp.parse(ds, row_meta_only=True))
    assert n_rid == n_cid, (
        "Input similarity matrix must be square.\n"
        "{} has dimensions {}x{}").format(ds, n_cid, n_rid)

    meta_df = get_col_meta(ds, external_metadata_id_field,
                           external_metadata_path, match_fields)

    # Add match_field column to meta_df after grouping by match_fields
    add_col_to_meta_df(meta_df, match_fields, MATCH_FIELD_NAME, SEPARATOR)

    sim_mat_result_record_list = create_sim_mat_result_records(
        ds, subset_fields, meta_df,
        external_metadata_path,
        external_metadata_id_field)

    for record in sim_mat_result_record_list:
        if len(sim_mat_result_record_list) > 1:
            logger.info("Evaluating sets, subsetted on: {}".format(record.subset_field_combo))
        record.results_df, record.nulls_gct = evaluate_sets_on_single_sim_mat(
            ds, set_definitions, aggregation_method, null_size, sets_per_chunk,
            low_memory_mode, max_set_size, record, meta_df)

    return sim_mat_result_record_list


def get_col_meta(ds, external_metadata_id_field, external_metadata_path,
                 match_fields):
    """ Get column metadata for input similarity matrix.
    Args:
        ds (string)
        external_metadata_id_field (string)
        external_metadata_path (string)
        match_fields (list of strings)

    Returns:
        meta_df (pandas df)
    """

    # Read in external metadata if provided
    if external_metadata_path is not None:
        meta_df_tmp = pd.read_csv(external_metadata_path, sep="\t")
        gct_col_meta = cpp.parse(ds, col_meta_only=True)

        # Make sure that each sample in GCT(X) is also in the external
        # metadata file, and that external_metadata_id_field uniquely
        # identifies samples
        check_external_meta(gct_col_meta, meta_df_tmp, external_metadata_id_field)

        # Set external_metadata_id_field as id
        meta_df_tmp.set_index(external_metadata_id_field, inplace=True)

        # Ignore metadata for things that are not in the dataset
        meta_df = meta_df_tmp.loc[gct_col_meta.index, :]

    else:
        meta_df = cpp.parse(ds, col_meta_only=True)
        assert not meta_df.empty, (
            "If external metadata not provided, column metadata must be "
            "embedded in the GCT(X) file.")

    return meta_df


def check_external_meta(gct_col_meta, external_meta, sample_id_field):
    """ Make sure that each sample in gct_col_meta is also in external_meta,
    and that sample_id_field uniquely identifies samples.

    Args:
        gct_col_meta (pandas df)
        external_meta (pandas df)
        sample_id_field (string)

    Returns:
        None

    """
    assert sample_id_field in external_meta.columns, (
        ("sample_id_field must be in external_meta.columns. sample_id_field: " +
         "{}, external_meta.columns: {}").format(
            sample_id_field, external_meta.columns.values))

    # Check that each sample in GCT(X) also in meta_df
    is_sample_in_external_bool_arr = gct_col_meta.index.isin(external_meta[sample_id_field])
    are_all_samples_in_external = all(is_sample_in_external_bool_arr)

    if not are_all_samples_in_external:
        first_sample_not_in_external_meta = gct_col_meta.index[~is_sample_in_external_bool_arr][0]
        errmsg = ("Not all samples in the GCT(X) are in the external " +
                  "metadata file. First such example: {}").format(first_sample_not_in_external_meta)
        raise AssertionError(errmsg)

    # Check also that sample_id_field identifies samples uniquely
    if not external_meta[sample_id_field].is_unique:
        duplicate_sample_ids = external_meta.loc[external_meta[sample_id_field].duplicated(), sample_id_field].values
        errmsg = ("sample_id_field {} must uniquely identify samples. " +
                  "Duplicate sample ids: {}").format(sample_id_field, duplicate_sample_ids)
        raise AssertionError(errmsg)


def add_col_to_meta_df(meta_df, match_fields, new_col_name, separator):
    """ Create a new column in the metadata df that will be used for
    identifying sets. Can use one or more match_fields to create this new
    column. If match_fields is empty, then the ids will be used.

    Args:
        meta_df (pandas df)
        match_fields (list of strings): fields to concatenate for the new column
        new_col_name (string): name of new column
        separator (string): separator used if concatenating fields

    Returns:
        meta_df (pandas df): modified in place with a new column added by the
            name of new_col_name

    """
    # Just use ids to create new column
    if not match_fields:
        meta_df[new_col_name] = meta_df.index.astype(str)

    # Otherwise, create a new column, possibly concatenating several fields
    else:
        for f in match_fields:
            assert f in meta_df.columns, (
                "match_field {} not in meta_df headers: {}".format(f, meta_df.columns))

        meta_df[new_col_name] = meta_df[match_fields].astype(str).apply(
            lambda x: x.str.cat(sep=SEPARATOR), axis=1)


def create_sim_mat_result_records(ds, subset_fields, meta_df,
                                  external_metadata_path,
                                  external_metadata_id_field):
    """ Create a list of similarity matrix result records, populated with
    one record for each subset combination of input similarity matrix. In
    the base case of no subsetting, this will be a list of one record.

    Args:
        ds (string)
        subset_fields (list of strings)
        meta_df (pandas df)
        external_metadata_id_field (string)

    Returns:
        result_record_list (list of SimMatResultRecord objects): each record
            will have the following attributes set:
                * sample_ids
                * subset_field_combo
                * subset_fname_suffix
    """
    result_record_list = []

    if len(subset_fields) > 0:
        subset_field_values_list = [list(meta_df[s_f_name].unique()) for s_f_name in subset_fields]
        subset_field_combos = list(itertools.product(*subset_field_values_list))
        gct_ids = meta_df.index
        for subset_combo in subset_field_combos:
            subset_ids = get_subset_ids(subset_combo, subset_fields,
                                        meta_df, external_metadata_id_field,
                                        gct_ids)
            if len(subset_ids) > 0:
                result_record_list.append(SimMatResultRecord(subset_ids, subset_combo))
    else:
        result_record_list.append(SimMatResultRecord(list(meta_df.index), ("")))

    return result_record_list


def get_subset_ids(subset_combo, subset_fields, meta_df,
                   external_metadata_id_field, gct_ids):
    """ Create a list of ids to subset from input similarity matrix, based on
    subset field combination criteria.
    For example, if subset_combo == ("A549", "24 h"), we will return a list of
    ids matching both those subset criteria.

    Args:
        subset_combo (tuple of strings)
        subset_fields (list of strings)
        meta_df (pandas df)
        external_metadata_id_field (string)
        gct_ids (list of strings)

    Returns:
        subset_ids (list of strings): list of ids to subset from similarity matrix
    """

    # we should only be here if we actually want to subset the sim mat
    assert len(subset_fields) > 0

    subset_masks = []
    for ii, subset_field_name in enumerate(subset_fields):
        subset_masks.append(meta_df[subset_field_name].isin([subset_combo[ii]]))
    all_subset_fields_mask = np.logical_and.reduce(subset_masks)

    # subsetted ids from metadata
    subset_ids = list(meta_df.loc[all_subset_fields_mask, :].index)
    # # ensure ids we want to subset on are in our input similarity matrix
    # subset_ids = [s_id for s_id in meta_subset_ids if s_id in gct_ids]

    return subset_ids


def evaluate_sets_on_single_sim_mat(ds, set_definitions, aggregation_method,
                                null_size, sets_per_chunk, low_memory_mode,
                                max_set_size, sim_mat_result_record, meta_df):
    """ Evaluate sets in a single similarity matrix. Return aggregated
    similarity, p-value, q-value, and n for each set in the similarity matrix.
    Also return the aggregated similarities for null distributions.

    Args:
        ds (string)
        set_definitions (string)
        aggregation_method (string)
        null_size (integer)
        sets_per_chunk (integer)
        low_memory_mode (bool)
        max_set_size (integer)
        sim_mat_result_record (SimMatResultRecord)
        meta_df (pandas df)

    Returns:
        results_df (pandas df): index is set names, columns are the following:
            agg_sim: aggregated similarity
            set_size: number of samples in this set
            pval: p-value for this set
            qval: q-value for this set
        nulls_gct (GCToo object): # rows = # of null iterations, # columns =
            # of unique real set sizes

    """

    # Subset metadata to sample ids
    meta_df = meta_df[meta_df.index.isin(sim_mat_result_record.sample_ids)]

    if low_memory_mode:
        if str.lower(os.path.splitext(ds)[1]) == ".gct":
            logger.warning(
                ("No point in running in low-memory mode if using a GCT, " +
                 "rather than GCTX, file. Will run in ordinary mode."))
            low_memory_mode = False
            gct = cpp.parse(ds, rid=sim_mat_result_record.sample_ids,
                        cid=sim_mat_result_record.sample_ids)
        else:
            gct = None

    else:
        gct = cpp.parse(ds, rid=sim_mat_result_record.sample_ids,
                    cid=sim_mat_result_record.sample_ids)

    # Read in GMT file if provided
    if set_definitions is not None:
        set_to_samples_dict_tmp = make_sets_from_gmt(set_definitions, meta_df, MATCH_FIELD_NAME)

    else:
        # Return a dictionary mapping sets to the samples comprising the sets
        set_to_samples_dict_tmp = make_sets_from_meta_df(meta_df, MATCH_FIELD_NAME)

    # Remove sets over a certain size
    set_to_samples_dict = remove_big_sets(set_to_samples_dict_tmp, max_set_size)

    # Get unique set sizes; need this for making nulls
    unique_set_sizes = set([len(set_samples) for set_samples in set_to_samples_dict.itervalues()])

    # Make nulls by picking random sets
    null_set_to_samples_dict = make_null_sets(
        meta_df, unique_set_sizes, null_size)
    logger.info("There are {} sets corresponding to {} unique set sizes.".format(
        len(set_to_samples_dict), len(null_set_to_samples_dict)))

    # If in low-memory mode, evaluate chunks of sets at a time
    if low_memory_mode:
        logger.info("Running in low-memory mode.")

        # Evaluate chunks of null sets
        set_size_to_agg_sim_dict = evaluate_chunks_of_null_sets(
            ds, null_set_to_samples_dict, aggregation_method, sets_per_chunk)

        # Evaluate chunks of real sets
        results_df = evaluate_chunks_of_real_sets(
            ds, set_to_samples_dict, set_size_to_agg_sim_dict,
            aggregation_method, sets_per_chunk)

    # Otherwise, evaluate all sets at once
    else:
        # TODO: evaluate_all_sets will be reorganized to act on gctoo (no dicts)
        results_df, set_size_to_agg_sim_dict = evaluate_all_sets(
            gct.data_df, set_to_samples_dict, null_set_to_samples_dict, aggregation_method)

    # Return null distributions as GCToo object
    nulls_gct = null_dict_to_gctoo(set_size_to_agg_sim_dict)

    return results_df, nulls_gct


def make_sets_from_meta_df(meta_df, match_field):
    """ Create sets using metadata from the GCT(X).

    Args:
        meta_df (pandas df)
        match_field (string): which field to use for creating sets

    Returns:
        set_to_samples_dict (dict): keys are set names, values
            are sample identifiers

    """

    # Get indices corresponding to each set
    dict_of_indices = meta_df.groupby(match_field).indices

    # Convert into dict where keys are set names, values are the sample identifiers
    set_to_samples_dict = {}
    for k, v in dict_of_indices.iteritems():

        if len(v) > 1:
            set_to_samples_dict[k] = set(meta_df.index[v])
        else:
            logger.debug("Set {} only has 1 element. Ignoring...".format(k))

    assert len(set_to_samples_dict) > 0, "No sets were created!"

    return set_to_samples_dict


def make_sets_from_gmt(path_to_gmt, meta_df, match_field):
    """ Get set definitions from the GMT file. Then use meta_df to create a
    mapping between set names and the ids of samples belonging to that set.
    Note that multiple samples can map to the same set member (i.e. an
    individual item from the GMT file), so we need to go through each
    set_member to get what samples it corresponds to.

    Args:
        path_to_gmt (string): path to GMT file
        meta_df (pandas df): metadata, either external or embedded
        match_field (string): metadata field in meta_df to use
            for matching to GMT set members

    Returns:
        set_to_samples_dict (dict): keys are sets, values are set (i.e. unique
            lists) of sample ids
    """
    assert os.path.exists(path_to_gmt)
    this_gmt = gmt.read(path_to_gmt)

    # Mess with gmt to make it a simple dictionary
    set_to_set_members_dict = {}
    for g in this_gmt:
        set_to_set_members_dict[g["head"]] = g["entry"]

    # Multiple samples can map to the same set member, so we need to go through
    # each set_member to get what samples it corresponds to
    set_to_samples_dict = {}
    for this_set, set_members in set_to_set_members_dict.items():

        all_samples_in_this_set = []
        for set_member in set_members:

            # If set_member not present in meta_df, skip it
            if not any(meta_df[match_field] == set_member):
                msg = ("set_member {} from set {} not present in meta_df. " +
                       "Skipping...").format(set_member, this_set)
                logger.warning(msg)
                continue

            # If set_member is present, add its samples to the growing list of
            # samples for this set
            samples_for_this_set_member = meta_df.index[meta_df[match_field] == set_member].tolist()
            all_samples_in_this_set.extend(samples_for_this_set_member)

        # Only include the set if it has more than one sample
        if len(all_samples_in_this_set) > 1:
            # Store list of samples as a set for faster access later
            set_to_samples_dict[this_set] = set(all_samples_in_this_set)
        else:
            msg = "set {} has fewer than 2 elements. Skipping...".format(this_set)
            logger.warning(msg)

    assert len(set_to_samples_dict) > 0, "No sets were created!"

    return set_to_samples_dict


def remove_big_sets(sets, max_set_size):
    """ Remove sets above max_set_size.

    Args:
        sets (dict): keys are set names, values are sets (i.e. unique lists)
            of sample ids
        max_set_size (int): maximum set size

    Returns:
        sets_clean (dict)

    """
    sets_clean = {}

    for k, v in sets.iteritems():
        if len(v) <= max_set_size:
            sets_clean[k] = v

    if len(sets_clean) != len(sets):
        num_removed = len(sets) - len(sets_clean)
        logger.info(("{} sets were removed because they had more than {} " +
                     "members.").format(num_removed, max_set_size))

    return sets_clean


def make_null_sets(meta_df, unique_set_sizes, null_size):
    """ Randomly pick samples from meta_df to create null sets. Create N=null_size
    sets for each unique set size.

    Args:
        meta_df (pandas df): metadata about samples
        unique_set_sizes (list of integers): the unique set sizes in these data
        null_size (integer): number of times to randomly sample

    Returns:
        null_set_to_samples_dict (dict): keys are unique_set_sizes, values are
            lists of sets (i.e. unique lists) of sample_ids

    """
    # Get sample identifiers
    samples = meta_df.index.values

    # For each unique_set_size, pick that many samples randomly and do it N=null_size times
    null_set_to_samples_dict = {}

    logger.info("Creating null sets...")
    for set_size in unique_set_sizes:
        this_list = [None] * null_size

        for ii in xrange(int(null_size)):
            these_random_samples = random.sample(samples, set_size)
            this_list[ii] = set(these_random_samples)

        null_set_to_samples_dict[set_size] = this_list

    return null_set_to_samples_dict


def evaluate_all_sets(data_df, real_sets, null_sets, aggregation_method):
    """ Evaluate all real and null sets. This is the function to use when
    the data matrix is not huge.

    Args:
        data_df (pandas df)
        real_sets (dict): keys are set names, values are sets (i.e. unique
            lists) of sample ids in data_df
        null_sets (dict): keys are fake set names, values are lists of sets
            (i.e. unique lists) of sample ids in data_df
        aggregation_method (string): how to aggregate the self sims

    Returns:
        results_df (pandas df): index is set names, columns are the following:
            agg_sim: aggregated similarity
            set_size: number of samples in a set
            pval: p-value for a set
            qval: q-value for a set
        set_size_to_agg_sim_dict (dict): keys are set sizes, values are
            arrays of aggregated similarities


    """
    # Initialize output
    results_df = pd.DataFrame(np.nan, index=real_sets.keys(),
                              columns=[AGG_SIM_COLUMN_NAME, SET_SIZE_COLUMN_NAME, PVAL_COLUMN_NAME])

    # Evaluate null sets (slowest part)
    logger.info("Evaluating null sets...")
    set_size_to_agg_sim_dict = evaluate_null_sets(data_df, null_sets, aggregation_method)

    # Evaluate each real set
    logger.info("Evaluating {} real sets...".format(len(real_sets)))
    for set_name, set_elements in real_sets.iteritems():

        (this_agg_sim, this_n, this_pval) = evaluate_one_set(
            data_df, set_elements, set_size_to_agg_sim_dict,
            aggregation_method)

        # Insert into output df
        results_df.loc[set_name, :] = [this_agg_sim, this_n, this_pval]

    # Compute q-values
    results_df[QVAL_COLUMN_NAME] = convert_to_qvals(results_df[PVAL_COLUMN_NAME].values)

    # Sort by the index
    results_df.sort_index(axis=0, inplace=True)

    return results_df, set_size_to_agg_sim_dict


def evaluate_one_set(data_df, set_elements, null_dict, aggregation_method):
    """ Compute similarity, n, and p-value for one set. Null similarities must
    be provided.

    Args:
        data_df (pandas df)
        set_elements (list): ids in data_df
        null_dict (dict): keys are set sizes, values are
            arrays of aggregated similarities
        aggregation_method (string): how to aggregate the self sims

    Returns:
        this_agg_sim (float)
        this_n (int)
        this_pval (float)

    """
    # Get set size
    this_n = int(len(set_elements))

    # Compute aggregated similarity
    this_agg_sim = get_agg_sim(data_df, set_elements, aggregation_method)

    # Compute p-value
    assert this_n in null_dict.keys()
    null_sims =  null_dict[this_n]
    this_pval = compute_pval(this_agg_sim, null_sims)

    return this_agg_sim, this_n, this_pval


def evaluate_null_sets(data_df, null_sets, aggregation_method):
    """ Return aggregated similarity for null sets in a dictionary, where keys
    are unique set sizes and values are arrays of aggregated sims.

    Args:
        data_df (pandas df)
        null_sets (dict): keys are set sizes, values are lists of sets (i.e.
            unique lists) containing sample ids
        aggregation_method (string): how to aggregated extracted similarities

    Returns:
        set_size_to_agg_sim_dict (dict): keys are set sizes, values are
            arrays of aggregated similarities

    """
    # Initialize output
    set_size_to_agg_sim_dict = {}

    for set_size, list_of_sets in null_sets.iteritems():
        these_sims = [get_agg_sim(data_df, s, aggregation_method) for s in list_of_sets]
        set_size_to_agg_sim_dict[set_size] = np.array(these_sims)

    return set_size_to_agg_sim_dict


def evaluate_chunks_of_real_sets(ds, real_sets, null_dict, aggregation_method, sets_per_chunk):
    """ Evaluate all real sets, chunks at a time. This is used in low-memory mode.

    Args:
        ds (string): path to GCT(X)
        real_sets (dict): keys are set names, values are sets (i.e. unique
            lists) of sample ids in GCT(X)
        null_dict (dict): keys are set sizes, values are
            arrays of aggregated similarities
        aggregation_method (string): how to aggregate the self sims
        sets_per_chunk (integer): number of sets to evaluate at a time in low-
            memory mode

    Returns:
        results_df (pandas df): index is set names, columns are the following:
            agg_sim: aggregated similarity
            set_size: number of samples in this set
            pval: p-value for this set
            qval: q-value for this set

    """
    # Create chunks of set_names and their corresponding elements
    chunks_of_set_names = chunkify(real_sets.keys(), sets_per_chunk)
    chunks_of_sets = chunkify(real_sets.values(), sets_per_chunk)
    chunks = zip(chunks_of_set_names, chunks_of_sets)

    list_of_dfs = []

    # Return aggregated similarities for sets in a single chunk
    # A chunk consists of tuples of set names and sets
    for chunk_num, chunk in enumerate(chunks):

        # Unpack
        set_names = chunk[0]
        list_of_sets = chunk[1]

        # Initialize output for this chunk
        this_df = pd.DataFrame(np.nan, index=set_names,
                               columns=[AGG_SIM_COLUMN_NAME, SET_SIZE_COLUMN_NAME, PVAL_COLUMN_NAME])

        logger.info("For real sets, evaluating chunk {} of {}.".format(
            chunk_num+1, len(chunks)))

        # Read in slice of GCT(X) corresponding to this chunk
        chunk_ids = get_chunk_ids(list_of_sets)
        g = cpp.parse(ds, rid=chunk_ids, cid=chunk_ids)

        for this_set_tuple in zip(*chunk):

            # Unpack tuple to get set name and set elements
            set_name = this_set_tuple[0]
            set_elements = this_set_tuple[1]

            (this_agg_sim, this_n, this_pval) = evaluate_one_set(
                g.data_df, set_elements, null_dict, aggregation_method)

            # Insert into output df
            this_df.loc[set_name, :] = [this_agg_sim, this_n, this_pval]

        # Append results for this chunk
        list_of_dfs.append(this_df)

    # Concatenate chunk results
    results_df = pd.concat(list_of_dfs, axis=0)

    # Compute q-values
    results_df[QVAL_COLUMN_NAME] = convert_to_qvals(results_df[PVAL_COLUMN_NAME].values)

    # Sort
    results_df.sort_index(axis=0, inplace=True)

    return results_df


def evaluate_chunks_of_null_sets(ds, null_sets, aggregation_method, sets_per_chunk):
    """ Return aggregated similarity for null sets in a dictionary, where keys
    are unique set sizes and values are arrays of aggregated sims. This is the
    version to use in low-memory mode.

    Args:
        ds (string): path to GCT(X)
        null_sets (dict): keys are set sizes, values are lists of sets (i.e.
            unique lists) containing sample ids
        aggregation_method (string): how to aggregated extracted similarities
        sets_per_chunk (integer): number of sets to evaluate at a time in low-
            memory mode

    Returns:
        set_size_to_agg_sim_dict (dict): keys are set sizes, values are
            arrays of aggregated similarities

    """
    # Initialize output
    set_size_to_agg_sim_dict = {}

    # Loop over set sizes
    for set_size, list_of_sets in null_sets.iteritems():

        # Initialize output for this set_size
        agg_sims_for_this_set_size = []

        # Create chunks of sets
        chunks_of_sets = chunkify(list_of_sets, sets_per_chunk)

        # Return aggregated similarities for sets in a single chunk
        # A chunk is a (smaller) list of sets
        for chunk_num, chunk in enumerate(chunks_of_sets):

            logger.info("For set size {}, evaluating chunk {} of {}.".format(
                set_size, chunk_num+1, len(chunks_of_sets)))

            # Read in slice of GCT(X) corresponding to this chunk
            chunk_ids = get_chunk_ids(chunk)
            g = cpp.parse(ds, rid=chunk_ids, cid=chunk_ids)

            # Get aggregated similarities
            these_sims = [get_agg_sim(g.data_df, s, aggregation_method) for s in chunk]
            agg_sims_for_this_set_size.extend(these_sims)

        assert len(agg_sims_for_this_set_size) == len(list_of_sets)
        set_size_to_agg_sim_dict[set_size] = np.array(agg_sims_for_this_set_size)

    return set_size_to_agg_sim_dict


def chunkify(my_list, items_per_chunk):
    """ Split a list up into chunks with n=items_per_chunk. If items_per_chunk
    is greater than the length of the list, return the list.

    Args:
        my_list (list)
        items_per_chunk (integer)

    Returns:
        chunks (list of lists)

    """
    if items_per_chunk > len(my_list):
        return [my_list]

    else:
        chunks = [my_list[ii:ii + items_per_chunk] for ii in xrange(
            0, len(my_list), items_per_chunk)]
        return chunks


def get_chunk_ids(list_of_sets):
    """ Extract all ids from list_of_sets.

    Args:
        list_of_sets (list): each set contains sample ids

    Returns:
        all_ids (list): all unique ids in list_of_sets

    """
    all_ids = []
    for s in list_of_sets:
        all_ids.extend(s)

    return list(set(all_ids))


def get_agg_sim(df, set_ids, aggregation_method):
    """ Return single value summarizing the self-similarities corresponding to
    these set_ids.

    Args:
        df (pandas df)
        set_ids (set of strings)
        aggregation_method (string)

    Returns:
        aggregated (float)

    """
    # Get rows and columns corresponding to set_ids
    small_df = df.loc[set_ids, set_ids]

    # Mask the lower triangle
    masked = small_df.where(np.triu(np.ones(small_df.shape), k=1).astype(np.bool))
    masked_vals = masked.values.flatten()

    # Make sure we have some non-NaN values
    if all(np.isnan(masked_vals)):
        aggregated = np.nan

    # Aggregate
    else:
        if aggregation_method == "q75":
            aggregated = np.nanpercentile(masked_vals, 75)

        elif aggregation_method == "median":
            aggregated = np.nanmedian(masked_vals)

        elif aggregation_method == "median_of_medians":
            small_df_matrix = small_df.values.astype("float")
            np.fill_diagonal(small_df_matrix, np.nan)
            aggregated = np.nanmedian(np.nanmedian(small_df_matrix, axis=1))

    return aggregated


def compute_pval(this_agg_sim, null_sims):
    """ Compute p-value by seeing how many similarities in the null are
     greater than the observed similarity.

    Args:
        this_agg_sim (float)
        null_sims (numpy array)

    Returns:
        pval

    """
    pval = np.divide((null_sims > this_agg_sim).sum(), float(len(null_sims)))
    return pval


def convert_to_qvals(pvals):
    """ Convert p-values to q-values using the Bonferroni-Hochberg approach.

    Args:
        pvals (numpy array)

    Returns:
        qvals (numpy array)

    """
    # Initialize output numpy array
    qvals = np.full(pvals.shape, np.nan)

    # Create mask to exclude missing values
    mask = np.isfinite(pvals)

    # Compute q-values
    qvals[mask] = multicomp.multipletests(pvals[mask], method="fdr_bh")[1]

    return qvals


def null_dict_to_gctoo(set_size_to_agg_sim_dict):
    """ Convert dictionary of null distributions to a GCToo object.

    Args:
        set_size_to_agg_sim_dict (dictionary): keys are set sizes; values
            are aggregated similarities corresponding to randomly sampled sets

    Returns:
        nulls_gct (GCToo object): empty row and column metadata dfs;
            columns of df are set sizes, rows are an integer index

    """
    # Make sure each entry in the dict has the same length; otherwise, pd.DataFrame will fail
    expected_length = len(set_size_to_agg_sim_dict[set_size_to_agg_sim_dict.keys()[0]])
    all_same_length = all([len(val) == expected_length for val in set_size_to_agg_sim_dict.itervalues()])
    assert all_same_length, (
        "All distributions in set_size_to_agg_sim_dict must have the same " +
        "length, but they don't. [len(val) for val in set_size_to_agg_sim_dict." +
        "itervalues()] : {}".format(
            [len(val) for val in set_size_to_agg_sim_dict.itervalues()]))

    # Create GCToo
    nulls_gct = GCToo.GCToo(pd.DataFrame(set_size_to_agg_sim_dict))

    return nulls_gct


def collate_results_dfs(sim_mat_result_record_list, subset_fields):
    """ Concatenate results dfs and annotate additional columns to specify
    how each result was subsetted.
    e.g. if subset_fields == ("cell_id", "pert_itime"), the resulting
    collated df might look like:
        agg_sim  set_size  pval  qval  cell_id  pert_itime
        ...                            A549     6 h
        ...                            ...
        ...                            A549     24 h
        ...                            ...
        ...                            MCF7     6 h
        ...                            ...
        ...                            MCF7     24 h
        ...                            ...      ...

    Args:
        sim_mat_result_record_list (list of SimMatResultRecord)
        subset_fields (list of strings)
    Returns:
        collated_results_df (pandas df): collated results dataframe to write to
            file, with subsetting parameters annotated as additional columns with
            subset fields as headers. No additional columns will appear if not
            subsetting.
    """
    if len(subset_fields) > 0:
        annotated_df_list = []
        for record in sim_mat_result_record_list:
            annotated_df = record.results_df.copy()
            for ii, subset_field_name in enumerate(subset_fields):
                annotated_df[subset_field_name] = record.subset_field_combo[ii]
            annotated_df_list.append(annotated_df)
        collated_results_df = pd.concat(annotated_df_list)
        collated_results_df.sort_values(subset_fields).sort_index()
    else:
        assert len(sim_mat_result_record_list) == 1, ("No fields to subset on, "
                                                      "but {} similarity matrix"
                                                      "result records were created. "
                                                      "Cannot collate multiple "
                                                      "results dfs when there are"
                                                      "no fields to subset on.")
        collated_results_df = sim_mat_result_record_list[0].results_df

    return collated_results_df


def write_null_gcts(sim_mat_result_record_list, out_dir):
    """ Save null distributions for each similarity matrix that eval_sets was
    run on, with size of matrix and any subset params as part of file name.

    Args:
        sim_mat_result_record_list (list of SimMatResultRecord)
        out_dir (string)

    Returns:
        nulls_fname_list (list of strings): list of written nulls file names
    """
    nulls_fname_list = []
    for record in sim_mat_result_record_list:
        full_out_nulls_name = os.path.join(
            out_dir,
            NULL_DIST_OUT_NAME.format(n_cols=record.nulls_gct.data_df.shape[1],
                                      n_rows=record.nulls_gct.data_df.shape[0],
                                      subset_fname_suffix=record.subset_fname_suffix))
        wgx.write(record.nulls_gct, full_out_nulls_name)
        nulls_fname_list.append(full_out_nulls_name)

    return nulls_fname_list


def write_yaml(args, full_out_text_name, nulls_fname_list):
    """ Save settings of evaluate_sets run, along with output filenames
    to a YAML file

    Args:
        args (dict)
        full_out_text_name (string)
        nulls_fname_list (list of strings)
    Returns:
        yaml_out_fname (string): name of written yaml file
    """
    yaml_out_fname = os.path.join(args.out_dir, YAML_FILE_OUT_NAME)
    with open(yaml_out_fname, "w") as yaml_handle:
        yaml_handle.write("results_fname: {}\n".format(full_out_text_name))
        yaml_handle.write("nulls_fnames: \n")
        for nulls_fname in nulls_fname_list:
            yaml_handle.write(" - {}\n".format(nulls_fname))
        for argname, argvalue in vars(args).iteritems():
            if type(argvalue) == list:
                yaml_handle.write("{}: \n".format(argname))
                for subvalue in argvalue:
                    yaml_handle.write(" - {}\n".format(subvalue))
            else:
                yaml_handle.write("{}: {}\n".format(argname, argvalue))
    return yaml_out_fname


if __name__ == "__main__":
    main()


