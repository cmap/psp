import sys
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import argparse
import logging

import cmapPy.pandasGEXpress.setup_GCToo_logger as setup_logger
import cmapPy.pandasGEXpress.parse as parse

logger = logging.getLogger(setup_logger.LOGGER_NAME)


def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required args
    parser.add_argument("--list_of_gcts", "-l", nargs="+", required=True,
                        help="space separated filepaths to 1+ input GCTs")

    parser.add_argument("--metadata_field", "-m", default="det_well_enrichment_score",
                        help="name of metadata field to plot on x axis")

    parser.add_argument("--output_name", "-o", default="probe_scatter.pdf",
                        help="name of output pdf file generated")

    parser.add_argument("-verbose", "-v", action="store_true", default=False,
                        help="increase the number of messages reported")

    return parser


def main(args):
    # Read GCTs into a list
    gctoo_list = [parse.parse(gct) for gct in args.list_of_gcts]

    # Create superset of all probes in GCTs
    probe_superset = create_probe_superset(gctoo_list)

    # Create pdf in which each page is a probe of the superset
    create_output_pdf(probe_superset, gctoo_list, args.metadata_field, args.output_name)


def create_probe_superset(gctoo_list):
    # Create list of sets of probes in each gct and return union of all sets
    list_of_probe_sets = [set(gct.data_df.index) for gct in gctoo_list]
    probe_superset = reduce(lambda a, b: a.union(b), list_of_probe_sets)

    return probe_superset


def create_output_pdf(probe_superset, gctoo_list, metadata_field, output_name):
    with PdfPages(output_name) as pdf:
        for probe in probe_superset:
            page_figure = plotify(probe, gctoo_list, metadata_field)
            pdf.savefig(page_figure)
            plt.close()

    return


def plotify(probe, gctoo_list, metadata_field):
    """ Iterates through provided GCTs to plot GCT values for given
    metadata field against probe quant values

    Args:
        probe (string) name of probe row
        gctoo_list (list of GCToo objects)
        metadata_field (string) name of metadata column

    Returns:
        figure (plot)
    """

    if len(gctoo_list) > 1:
        plt.figure()
        fig, axes = plt.subplots(1, len(gctoo_list), sharey=True, sharex=True)
        plt.suptitle(probe, fontsize=16)
        plt.xlabel(metadata_field)

        for i in range(len(gctoo_list)):
            gct = gctoo_list[i]
            x_vals = gct.col_metadata_df.loc[:, metadata_field]

            # Account for GCTs in which probe field may have been filtered
            try:
                y_vals = gct.data_df.loc[probe, :]
            except KeyError as error:
                # If probe does not exist in GCT y values are null
                y_vals = pd.Series(index=gct.data_df.columns)

            # Set up plot sizing
            axes[i].tick_params(axis='both', which='major', labelsize=8)
            axes[i].set_title(gct.src, fontsize=5)

            axes[i].scatter(x_vals, y_vals)

            # Set y axis label on first / left-most plot only
            if i == 0:
                axes[i].set_ylabel("Probe Quant Value")
    else:

        fig = plt.figure()
        plt.title(probe)
        plt.xlabel(metadata_field)
        plt.ylabel("Probe Quant Value")
        x_vals = gctoo_list[0].col_metadata_df.loc[:, metadata_field]

        try:
            y_vals = gctoo_list[0].data_df.loc[probe, :]
        except KeyError as error:
            # If probe does not exist in GCT y values are null
            y_vals = pd.Series(index=gctoo_list[0].data_df.columns)

        plt.scatter(x_vals, y_vals)

    plt.close()
    return fig


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    setup_logger.setup(verbose=args.verbose)
    logger.debug("args: {}".format(args))

    main(args)