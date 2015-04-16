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


def bismark_analysis(input_cov_list, min_cov=10, fasta_ref=None, plot_scatters=None):
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

    # Build nicer sample names
    sample_names = {}
    for fn in input_cov_list:
        name = os.path.basename(fn)
        remove = ['_val_1', '_val_2', '.fq', '.gz', '_bismark', '_bt2', '_bt1', '_pe', '_se', '.deduplicated', '.bismark', '.gwCov', '.cov']
        for r in remove:
            name = name.replace(r, '')
        sample_names[fn] = name

    # Plot dendrogram
    make_dendrogram(data, sample_names)

    # Auto - only plot scatters if 4 or less samples (24 plots)
    if plot_scatters is None:
        samps = len(input_cov_list)
        num_scatters = (samps**2)-((samps**2)/2)-(samps/2)
        if num_scatters <= 24:
            plot_scatters = True
            logging.info("{} samples given, plotting {} scatter plots".format(samps, num_scatters))
        else:
            plot_scatters = False
            logging.info("{} samples would give {} scatter plots, so skipping. Use '--scatter t' to force.".format(samps, num_scatters))

    # Plot scatter plots
    if plot_scatters is True:
        for i in range(0, len(input_cov_list)-1):
            for j in range(i+1, len(input_cov_list)):
                x = input_cov_list[j]
                y = input_cov_list[i]
                plot_meth_scatter(data[x], data[y], sample_names[x], sample_names[y])



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



def make_dendrogram(data, sample_names=False, output_fn=False, plot_title=False):
    """
    Plots a hierarchical clustering dendrogram
    """

    # Output filename
    if output_fn is False:
        output_fn = "sample_dendrogram"

    # Find shared keys present in all datasets
    shared_keys = list(set.intersection(* map(set, data.values())))

    # How many did we see in total?
    all_keys = list(set.union(* map(set, data.values())))

    # Count and log - vars used later in plot too
    num_total = len(all_keys)
    num_shared = len(shared_keys)
    skipped_keys = num_total - num_shared
    percent_shared = (float(num_shared)/float(num_total))*100
    percent_skipped = (float(skipped_keys)/float(num_total))*100
    logging.info("Found {:.2e} ({:.0f}%) cytosines present in all datasets for dendrogram. Skipped {:.2e} ({:.0f}%) out of {:.2e} total observed.".format(num_shared, percent_shared, skipped_keys, percent_skipped, num_total))

    # Build data array
    df = []
    dslabels = []
    for fn, d in data.iteritems():
        df.append([d[k]['methylation'] for k in shared_keys])
        if sample_names is not False:
            dslabels.append(sample_names[fn])
        else:
            dslabels.append(fn)

    # Calculate linkage and plot dendrogram
    fig = plt.figure(figsize=(12,4.5), tight_layout={'rect':(0,0.04,1,1)})
    axes = plt.axes(frameon=False)
    link = cluster.hierarchy.linkage(df, metric='correlation')
    # Above threshold colour is hardcoded as blue. Valls put in a PR to scipy for this.
    # Dev version has the attr now, let's use it if we can..
    try:
        den = cluster.hierarchy.dendrogram(link, count_sort=True, labels=dslabels, above_threshold_color='#AAAAAA')
    except TypeError:
        den = cluster.hierarchy.dendrogram(link, count_sort=True, labels=dslabels)

    # Labels
    if plot_title is False:
        plot_title = "Sample Clustering by Methylation"
    plt.title(plot_title)

    # Formatting
    axes.tick_params(right=False, axis='both', labelsize=6)
    for label in axes.get_xticklabels():
        label.set_rotation(90)

    # Add a label about number of Cs used
    plt.figtext(0.99, 0.03, r'Calculated using {:.2e} cytosines found in all samples ({:.0f}% of all seen)'.format(num_shared, percent_shared),
                horizontalalignment='right', fontsize=6, transform=axes.transAxes)

    # Save Plots
    png_fn = "{}.png".format(output_fn)
    pdf_fn = "{}.pdf".format(output_fn)
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)

    # Return the filenames
    return {'png': png_fn, 'pdf': pdf_fn}



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

    # Set up variables
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
    fig = plt.figure(figsize=(8,7))
    axes = fig.add_subplot(111, aspect=1.0)
    plt.hist2d(x_vals, y_vals, bins=200, norm=mpl.colors.LogNorm(), cmap=cm.Blues)
    cbar = plt.colorbar()

    # Axis labels and titles
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.title("CpG Methylation Percentages")
    cbar.set_label('Cytosine Counts')

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
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)

    # Return the filenames
    return {'png': png_fn, 'pdf': pdf_fn}



if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("Compare bismark methylation coverage data between samples")
    parser.add_argument("-f", "--fasta", dest="fasta_ref",
                        help="Genome Fasta reference for coverage analysis")
    parser.add_argument("-s", "--scatter", dest="plot_scatters", default='auto', choices=['auto', 't', 'f'],
                        help="Plot pairwise scatter plots? Auto = only if 4 or less samples (24 plots)")
    parser.add_argument("-l", "--log", dest="log_level", default='info', choices=['debug', 'info', 'warning'],
                        help="Level of log messages to display")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("input_cov_list", metavar='<cov file>', nargs="+",
                        help="List of input genome-wide coverage filenames")
    kwargs = vars(parser.parse_args())

    # Scatter var munging
    if kwargs['plot_scatters'] == 'auto': kwargs['plot_scatters'] = None
    if kwargs['plot_scatters'] == 't': kwargs['plot_scatters'] = True
    if kwargs['plot_scatters'] == 'f': kwargs['plot_scatters'] = False

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
