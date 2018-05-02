"""
plot_conns_for_one_target.py

This fcn should work for any connectivity metric. However, remember to
adjust OUT_SUFFIX. Also, if it's a metric where small is good (e.g. p-values),
you might want to switch the sort order to ascending=True.

"""

import pandas as pd
import os
import sys
import argparse
import parse_gct as pg
import matplotlib.pyplot as plt

# TODO(lev): use setup_logger

OUT_SUFFIX = "_sorted_conn.png"

def build_parser():
    """Build argument parser."""

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("in_gct_path", type=str, help="full path to input gct")
    parser.add_argument("out_dir", type=str, help="where to save output")
    parser.add_argument("out_prefix", type=str, help="prefix for naming output figure and its title")
    parser.add_argument("target_id", type=str, help="which row of connectivity matrix to extract")

    parser.add_argument("-queries_to_highlight", "-qth", nargs="*", type=str, default=None,
                        help="which queries to highlight")
    parser.add_argument("-conn_metric", type=str, default="KS test statistic",
                        help="connectivity metric to use for plot labeling")

    return parser


def main(args):

    # Parse input gct
    print("gct_path: {}".format(args.in_gct_path))
    in_gct = parse.parse(args.in_gct_path)

    # Extract target_id row
    target_s = in_gct.data_df.loc[args.target_id, :]
    assert isinstance(target_s, pd.Series)

    # Sort conns
    target_s_sorted = target_s.sort_values(ascending=True)

    # Return the number of queries (i.e. number of points on plot)
    n_queries = target_s_sorted.shape[0]
    print("number of queries: {}".format(n_queries))

    # Create plot
    fig1 = plt.figure()
    plt.hold(True)
    plotting_inds = range(1, n_queries + 1)
    plt.scatter(target_s_sorted, plotting_inds)

    # Add labels and title
    plt.ylabel('Query')
    plt.xlabel(args.conn_metric)
    plt.title((args.out_prefix + ", target_id: {}").format(args.target_id))
    plt.ylim([0, n_queries + 1])
    plt.setp(plt.gca(), "yticklabels", [])

    # Replot and add labels for just queries_to_highlight
    if args.queries_to_highlight is not None:
        n_queries_to_highlight = len(args.queries_to_highlight)

        for query_id in args.queries_to_highlight:
            query_ind = target_s_sorted.index.get_loc(query_id)
            conn_value = target_s_sorted[target_s_sorted.index == query_id][0]
            plt.plot(conn_value, query_ind, "or")
            plt.annotate(query_id, (conn_value, query_ind),
                         horizontalalignment="left",
                         verticalalignment="top")

    # Save figure
    save_name = os.path.join(args.out_dir, str.lower((args.out_prefix + OUT_SUFFIX)))
    plt.savefig(save_name)


if __name__ == "__main__":
    args = build_parser().parse_args(sys.argv[1:])
    main(args)

