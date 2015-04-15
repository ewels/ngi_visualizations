#!/usr/bin/python
"""
count_biotypes.py
Takes aligned reads and uses HTSeq to count overlaps with different
biotype feature flags. Plots biotype proportions and a histogram
of read lengths broken up by biotype.
Input: GTF annotation file and aligned BAM files
"""

from __future__ import print_function

import argparse
from collections import defaultdict
import logging
import numpy as np
import os
from scipy import cluster, stats

# Import matplot lib but avoid default X environment
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm # colours


def bismark_analysis(input_cov_list, min_cov=10):
    """
    Run bismark comparison analysis - main function
    """
    # Load the data
    data = {}
    for fn in input_cov_list:
        data[fn] = load_bismark_gwCov(fn, min_cov)
        if data is None:
            logging.error("Could not parse coverage file: {} Skipping..".format(fn))
            del (data[fn])

    # Plot dendrogram
    make_dendrogram(data)
    os.exit(0)

    # Plot scatter plots
    for i in range(0, len(input_cov_list)-1):
        for j in range(i+1, len(input_cov_list)):
            x = input_cov_list[j]
            y = input_cov_list[i]
            x_name = os.path.basename(x)
            y_name = os.path.basename(y)
            remove = ['_val_1', '_val_2', '.fq', '.gz', '_bismark', '_bt2', '_bt1', '_pe', '_se', '.deduplicated', '.bismark', '.gwCov', '.cov']
            for r in remove:
                x_name = x_name.replace(r, '')
                y_name = y_name.replace(r, '')
            plot_meth_scatter(data[x], data[y], x_name, y_name)



def load_bismark_gwCov(fn, min_cov=10):
    """
    Load a Bismark coverage file into memory. Takes normal coverage
    or genome-wide coverage files. The former are much, *much* faster
    to load.
    """
    logging.info("Loading {}".format(fn))
    data = defaultdict(dict)
    i = 0
    total_cs = 56490324
    try:
        with open(fn, 'r') as fh:
            for line in fh:
                i += 1
                if i % 10000000 == 0:
                    logging.info("  ..line {0}M of approximately {1:.1f}M in human genome ({2:.2f}%)".format(i/1000000, float(total_cs)/1000000, (float(i)/float(total_cs))*100))

                # line = line.strip() # commented out for speedup
                c = line.split()

                # Genome wide coverage input
                if len(c) == 7:
                    meth = float(c[3])
                    unmeth = float(c[4])
                    cov = meth + unmeth
                    # strand = c[2]
                    methp = (meth / cov)*100 if (cov > 0) else None
                elif len(c) == 6:
                    meth = float(c[4])
                    unmeth = float(c[5])
                    cov = meth + unmeth
                    methp = c[3]
                else:
                    logging.error("Error: Coverage file {} had {} columns instead of 6 or 7.".format(fn, len(c)))
                    return None

                if cov >= min_cov:
                    # Add the coverage and % methylation scores
                    key = "{}_{}".format(c[0], c[1]) # chr_pos
                    data[key]['coverage'] = cov
                    data[key]['methylation'] = methp

    except IOError as e:
        logging.error("Error loading coverage file: {}".format(fn))
        raise IOError(e)

    return data



def plot_meth_scatter (sample_1, sample_2, x_lab, y_lab, output_fn=False):
    """
    Plots scatter plot of methylation score counts. Calculates r-squared
    correlation and spearman's correlation coefficient scores
    and prints these to STDOUT as well as including it in the graph.
    Returns a dict containing filenames of PNG and PDF graphs as a dict for
    each sample.
    """

    # Output filename
    if output_fn is False:
        output_fn = "{}_{}".format(y_lab, x_lab)

    # SET UP PLOT
    fig = plt.figure()
    axes = fig.add_subplot(111, aspect=1.0)
    x_vals = []
    y_vals = []
    all_x_vals = []
    all_y_vals = []
    missing_positions = 0
    none_positions = 0

    # Find intersection between samples
    shared_keys = sample_1.viewkeys() & sample_2.viewkeys()
    s1_shared = (float(len(shared_keys))/float(len(sample_1)))*100
    s2_shared = (float(len(shared_keys))/float(len(sample_1)))*100
    if len(shared_keys) == 0:
        logging.warn("Warning: No shared keys found for {} and {}. Skipping scatter plot.".format(x_lab, y_lab))
        return None
    logging.info("Found {} shared keys between {} and {} ({:.2f}% and {:.2f}%)".format(len(shared_keys), x_lab, y_lab, s1_shared, s2_shared))

    # Collect the methylation scores
    zero_percent = 0
    hundred_percent = 0
    for pos in shared_keys:
        try:
            thisy = float(sample_2[pos]['methylation'])
            thisx = float(sample_1[pos]['methylation'])
            all_y_vals.append(thisy)
            all_x_vals.append(thisx)
            if thisy == 100 and thisx == 100:
                hundred_percent +=1
            elif thisy == 0 and thisx == 0:
                zero_percent +=1
            else:
                y_vals.append(thisy)
                x_vals.append(thisx)
        except TypeError:
            none_positions += 1
    if missing_positions > 0:
        missing_p = (missing_positions/float(len(sample_1)))*100
        logging.warn("Warning: {} positions out of {} mentioned in {} not found in {} ({:.2f}%)".format(missing_positions, len(sample_1), x_lab, y_lab, missing_p))
    if none_positions > 0:
        logging.warn("Warning: {} positions skipped due to zero coverage in at least one sample for {} and {}".format(none_positions, x_lab, y_lab))
    if len(x_vals) == 0 or len(x_vals) != len(y_vals):
        logging.warn("Error: Bad lengths of lists for plotting(x = {}, y = {}) for {} and {}. Skipping.".format(len(x_vals), len(y_vals), x_lab, y_lab))
        return None

    # Calculate the r squared
    corr = np.corrcoef(all_x_vals, all_y_vals)[0,1]
    r_squared = corr ** 2
    print ("R squared for {} = {:.2f}".format(output_fn, r_squared))

    # Caculate the spearman's rank correlation coefficient
    rho, pvalue = stats.spearmanr(all_x_vals, all_y_vals)
    print ("Spearman's rank correlation coefficient for {} = {:.2f} (p < {:.2e})".format(output_fn, rho, pvalue))

    # Plot the histogram
    plt.hist2d(x_vals, y_vals, bins=200, norm=mpl.colors.LogNorm(), cmap=cm.Blues)
    cbar = plt.colorbar()

    # Axis labels and titles
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.title("CpG Methylation Percentages")
    cbar.set_label('Log Cytosine Counts')

    # Axes scales
    axes.set(xlim=[0,100])
    axes.set(ylim=[0,100])

    # Tidy axes
    axes.set_axisbelow(True)
    axes.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)
    cbar.ax.tick_params(axis='y', direction='out')

    # Add a label about r squared and spearman's
    plt.subplots_adjust(bottom=0.17)
    plt.text(1, -0.14, r'$r^2$ = {:.2f}'.format(r_squared),
                horizontalalignment='right', fontsize=8, transform=axes.transAxes)
    plt.text(1, -0.17, r'$\rho$ = {:.2f}'.format(rho),
                horizontalalignment='right', fontsize=8, transform=axes.transAxes)

    # Add a label about 0% and 100%
    plt.text(0, -0.14, r'Both 100% Methylated = {:.2e}'.format(zero_percent),
                horizontalalignment='left', fontsize=8, transform=axes.transAxes)
    plt.text(0, -0.17, r'Both 0% Methylated = {:.2e}'.format(hundred_percent),
                horizontalalignment='left', fontsize=8, transform=axes.transAxes)

    # SAVE OUTPUT
    png_fn = "{}.png".format(output_fn)
    pdf_fn = "{}.pdf".format(output_fn)
    logging.info("Saving to {}".format(png_fn))
    plt.savefig(png_fn)
    logging.info("Saving to {}".format(pdf_fn))
    plt.savefig(pdf_fn)

    # Return the filenames
    return {'png': png_fn, 'pdf': pdf_fn}



def make_dendrogram(data, output_fn=False):
    """
    Plots a hierarchical clustering dendrogram
    """

    # Output filename
    if output_fn is False:
        output_fn = "sample_dendrogram"

    # Count occurance of positions in datasets
    keys = {}
    for d in data.values():
        for k in d.keys():
            if k not in keys:
                keys[k] = 1
            else:
                keys[k] += 1

    shared_keys = 0
    num_discarded = 0
    for k,c in keys.iteritems():
        if c == len(data):
            shared_keys += 1
            # TODO - make list of data points here or something
        else:
            num_discarded += 1
    print ("{} shared, {} ignored".format(len(shared_keys), num_discarded))

    # Calculate linkage
    # linkage = cluster.hierarchy.linkage



if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("Compare bismark methylation coverage data between samples")
    parser.add_argument("-l", "--log", dest="log_level", default='info', choices=['debug', 'info', 'warning'],
                        help="Level of log messages to display")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("input_cov_list", metavar='<cov file>', nargs="+",
                        help="List of input genome-wide coverage filenames")
    kwargs = vars(parser.parse_args())

    # Initialise logger
    numeric_log_level = getattr(logging, kwargs['log_level'].upper(), None)
    if kwargs['log_output'] != 'stdout':
        logging.basicConfig(filename=kwargs['log_output'], format='', level=numeric_log_level)
    else:
        logging.basicConfig(format='', level=numeric_log_level)
    # Remove logging parameters
    kwargs.pop('log_level', None)
    kwargs.pop('log_output', None)

    # Call bismark_analysis()
    bismark_analysis(**kwargs)
