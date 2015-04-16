#!/usr/bin/python
"""
bismark.py
Takes Bismark coverage BED files (generated using bismark2bedGraph)
and creates a number of plots:
 * Clustering dendrogram of samples
 * Pairwise scatter plots
 * Coverage analysis (if given a fasta reference with --fasta)
Can be run on the command line or imported and used directly.
"""

from __future__ import print_function

import argparse
from collections import defaultdict
import logging
import numpy as np
import os
from scipy import cluster, stats
from Bio import SeqIO

# Import matplot lib but avoid default X environment
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm # colours


def bismark_analysis(input_cov_list, min_cov=10, fasta_fn=None, regions_fn=False, plot_scatters=None):
    """
    Run bismark comparison analysis - main function
    """
    # Load the data
    data = {}
    for fn in input_cov_list:
        data[fn] = load_bismark_cov(fn, min_cov)
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
    make_dendrogram(data, min_cov, sample_names)

    # Auto - only plot scatters if 4 or less samples (24 plots)
    if plot_scatters is None:
        samps = len(input_cov_list)
        num_scatters = (samps**2)-((samps**2)/2)-(samps/2)
        if num_scatters <= 10:
            plot_scatters = True
            logging.info("{} samples given, plotting {} scatter plots".format(samps, num_scatters))
        else:
            plot_scatters = False
            logging.info("{} samples would give {} scatter plots, so skipping. Use '--scatter t' to force.".format(samps, num_scatters))

    # Plot scatter plots and work out correlation scores
    rsquared = defaultdict(lambda: defaultdict(float))
    rho = defaultdict(lambda: defaultdict(float))
    for i in range(0, len(input_cov_list)-1):
        for j in range(i+1, len(input_cov_list)):
            x = input_cov_list[j]
            y = input_cov_list[i]
            results = plot_meth_scatter(data[x], data[y], sample_names[x], sample_names[y], min_cov, plot_scatters)
            rsquared[sample_names[x]][sample_names[y]] = results['r_squared']
            rho[sample_names[x]][sample_names[y]] = results['rho']

    # Do coverage analysis
    if fasta_fn is not None:
        if regions_fn is not None:
            regions = load_capture_regions(regions_fn)
        else:
            regions = None
        (fasta_ref, captured_cgs) = load_fasta_cpg(fasta_fn, regions)
        total_cg_count = len(fasta_ref)
        coverage_decay_plot(data, sample_names, total_cg_count, captured_cgs)


def load_bismark_cov(fn, min_cov=10):
    """
    Load a Bismark coverage file into memory. Takes normal coverage
    or genome-wide coverage files. The former are much, *much* faster
    to load.
    """
    logging.info("Loading {}".format(fn))
    data = defaultdict(dict)
    try:
        with open(fn, 'r') as fh:
            for line in fh:
                c = line.split()

                # Genome wide coverage input
                if len(c) == 7:
                    meth = float(c[3])
                    unmeth = float(c[4])
                    cov = meth + unmeth
                    # strand = c[2]
                    if c[5] is not 'CG':
                        next
                    methp = (meth / cov)*100 if (cov > 0) else None
                elif len(c) == 6:
                    meth = float(c[4])
                    unmeth = float(c[5])
                    cov = meth + unmeth
                    methp = c[3]
                else:
                    logging.error("Error: Coverage file {} had {} columns instead of 6 or 7.".format(fn, len(c)))
                    return None

                # Add the coverage and % methylation scores
                chrom = c[0]
                chrom = chrom.replace('chr','')
                key = "{}_{}".format(chrom, c[1]) # chr_pos
                data[key]['coverage'] = cov
                data[key]['methylation'] = methp

    except IOError as e:
        logging.error("Error loading coverage file: {}".format(fn))
        raise IOError(e)

    return data



def make_dendrogram(data, min_cov=10, sample_names=False, output_fn=False, plot_title=False):
    """
    Plots a hierarchical clustering dendrogram based on sample methylation
    scores (for positions seen in all samples)
    """

    # Output filename
    if output_fn is False:
        output_fn = "sample_dendrogram"

    # Remove any positions with less than the threshold observations
    for f in data.keys():
        for p in data[f].keys():
            if data[f][p]['coverage'] < min_cov:
                del(data[f][p])

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
    logging.info("Found {:.2e} ({:.0f}%) cytosines present in all datasets for dendrogram. Skipped {:.2e} ({:.0f}%) out of {:.2e} total observed above coverage threshold of {}.".format(num_shared, percent_shared, skipped_keys, percent_skipped, num_total, min_cov))

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
    plt.figtext(0.99, 0.03, r'Calculated using {:.2e} cytosines found in all samples ({:.0f}% of all seen with coverage > {})'.format(num_shared, percent_shared, min_cov),
                horizontalalignment='right', fontsize=6, transform=axes.transAxes)

    # Save Plots
    png_fn = "{}.png".format(output_fn)
    pdf_fn = "{}.pdf".format(output_fn)
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)

    # Return the filenames
    return {'png': png_fn, 'pdf': pdf_fn}



def plot_meth_scatter (sample_1, sample_2, x_lab, y_lab, min_cov=10, plot=True, output_fn=None):
    """
    Plots scatter plot of methylation score counts. Calculates r-squared
    correlation and spearman's correlation coefficient scores
    and prints these to STDOUT as well as including it in the graph.
    Returns a dict containing filenames of PNG and PDF graphs as a dict for
    each sample.
    """

    # Remove any positions with less than the threshold observations
    for p in sample_1:
        if sample_1[p]['coverage'] < min_cov:
            del(sample_1[p])
    for p in sample_2:
        if sample_2[p]['coverage'] < min_cov:
            del(sample_2[p])

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
    logging.debug("Found {} shared keys between {} and {} ({:.2f}% and {:.2f}%)".format(len(shared_keys), x_lab, y_lab, s1_shared, s2_shared))

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

    # Caculate the spearman's rank correlation coefficient
    rho, pvalue = stats.spearmanr(all_x_vals, all_y_vals)

    # Plot the histogram
    if plot:
        fig = plt.figure(figsize=(8,7))
        axes = fig.add_subplot(111, aspect=1.0)
        plt.hist2d(x_vals, y_vals, bins=200, norm=mpl.colors.LogNorm(), cmap=cm.Blues)
        cbar = plt.colorbar()

        # Axis labels and titles
        plt.xlabel(x_lab)
        plt.ylabel(y_lab)
        plt.title("CpG Methylation Percentages (Coverage >= {})".format(min_cov))
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
        if output_fn is None:
            output_fn = "{}_{}_methylation_scatter".format(y_lab, x_lab)
        png_fn = "{}.png".format(output_fn)
        pdf_fn = "{}.pdf".format(output_fn)
        plt.savefig(png_fn)
        plt.savefig(pdf_fn)

    # Not plotting
    else:
        png_fn = None
        pdf_fn = None

    # Return the filenames and correlation scores
    return {'png': png_fn, 'pdf': pdf_fn, 'r_squared': r_squared, 'rho': rho, 'rho_pvalue': pvalue }



def load_capture_regions(regions_fn):
    """
    Load a Fasta reference file and finds positions of all CpGs
    in the genome. Returned as a list of string in format chr_pos
    """

    logging.info("Loading capture regions reference {}".format(regions_fn))
    regions = defaultdict(lambda: defaultdict(int))
    try:
        with open(regions_fn, 'r') as fh:
            for line in fh:
                c = line.split()
                if len(c) != 4:
                    continue
                chrom = c[0]
                chrom = chrom.replace('chr','')
                regions[chrom][int(c[1])] = int(c[2])

    except IOError as e:
        logging.error("Error loading capture regions file: {}".format(regions_fn))
        raise IOError(e)

    return regions


def load_fasta_cpg(fasta_fn, cap_regions=None):
    """
    Load a Fasta reference file and finds positions of all CpGs
    in the genome. Returned as a list of string in format chr_pos
    """

    logging.info("Loading Fasta reference {}".format(fasta_fn))
    fasta_ref = list()

    cap_starts = defaultdict(str)
    if cap_regions is not None:
        captured_cgs = list()
        for chrom in cap_regions:
            cap_starts[chrom] = cap_regions[chrom].keys()
    else:
        captured_cgs = None

    fasta_sequences = SeqIO.parse(open(fasta_fn),'fasta')
    for fasta in fasta_sequences:
        chrom, sequence = fasta.id, str(fasta.seq)
        chrom = chrom.replace('chr','')
        if cap_regions is None:
            logging.info("  ..searching chromosome {} for CpGs - found {} so far".format(chrom, len(fasta_ref)))
        else:
            logging.info("  ..searching chromosome {} for CpGs - found {} so far, {} captured".format(chrom, len(fasta_ref), len(captured_cgs)))
        pos = sequence.find("CG")
        while pos >= 0:
            fasta_ref.append('{}_{}'.format(chrom, pos))
            pos = sequence.find("CG", pos + 1)
            if cap_regions is not None:
                if chrom in cap_starts:
                    # Find closest capture region to this CG position
                    start = min(cap_starts[chrom], key=lambda x:abs(x-pos))
                    if start > pos:
                        sindex = cap_starts[chrom].index(start)
                        start = cap_starts[chrom][sindex - 1]
                    end = cap_regions[chrom][start]

                    # Does our position overlap this region?
                    if pos >= start and pos <= end:
                        captured_cgs.append('{}_{}'.format(chrom, pos))

    return (fasta_ref, captured_cgs)


def coverage_decay_plot(data, sample_names=None, total_cg_count=False, captured_cgs=None, output_fn=None):
    """
    Plots a coverage decay plot, y axis as percentage of genome-wide CpGs
    and y axis increasing coverage thresholds
    """
    # Collect counts of each coverage
    for fn, d in data.iteritems():

        logging.info("Saving coverage display plot for {}".format(fn))

        # Count the occurance of both types of coverage
        coverages = defaultdict(int)
        captured_coverages = defaultdict(int)
        for pos, arr in d.iteritems():
            coverages[arr['coverage']] += 1
            if captured_cgs is not None and pos in captured_cgs:
                captured_coverages[arr['coverage']] += 1

        # Build the genome wide coverage plotting lists
        x_vals = []
        y_vals = []
        count = 0
        for c in range(int(max(coverages.keys())), 0, -1):
            count += coverages[c]
            c_percent = (float(count) / float(total_cg_count))*100
            x_vals.append(c)
            y_vals.append(c_percent)
        x_vals.append(0)
        y_vals.append(100)

        # Build the captured coverage plotting lists
        if captured_cgs is not None:
            cap_x_vals = []
            cap_y_vals = []
            count = 0
            total_captured_cgs = len(captured_cgs)
            for c in range(int(max(captured_coverages.keys())), 0, -1):
                count += captured_coverages[c]
                c_percent = (float(count) / float(total_captured_cgs))*100
                cap_x_vals.append(c)
                cap_y_vals.append(c_percent)
            cap_x_vals.append(0)
            cap_y_vals.append(100)

        # Draw the plot
        fig = plt.figure(figsize=(10,6))
        plt.plot(x_vals, y_vals, 'b-')
        if captured_cgs is not None:
            plt.plot(cap_x_vals, cap_y_vals, 'g--')

        # Tidy axes
        plt.xlim(0,50)
        plt.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)

        # Labels
        try:
            sample_name = sample_names[fn]
        except KeyError:
            sample_name = fn

        plt.xlabel('CpG Coverage')
        plt.ylabel('Percentage of CGs')
        plt.title("Coverage Decay Plot: {}".format(sample_name))

        # Save the plot
        if output_fn is None:
            output_fn = "{}_coverage_decay".format(sample_name)
        plt.savefig('{}.png'.format(output_fn))
        plt.savefig('{}.pdf'.format(output_fn))



if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("Compare bismark methylation coverage data between samples")
    parser.add_argument("-m", "--min-cov", dest="min_cov", default=10,
                        help="Minimum coverage to use for scatter plots and dendrogram")
    parser.add_argument("-f", "--fasta", dest="fasta_fn",
                        help="Genome Fasta reference for coverage analysis")
    parser.add_argument("-c", "--capture", dest="regions_fn",
                        help="BED file containing capture regions")
    parser.add_argument("-s", "--scatter", dest="plot_scatters", default='auto', choices=['auto', 't', 'f'],
                        help="Plot pairwise scatter plots? Auto = only if 10 or less plots")
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
