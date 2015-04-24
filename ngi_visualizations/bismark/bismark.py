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

# Set the default color cycle
mpl.rcParams['axes.color_cycle'] = ['#AAAAAA', '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#fb9a99', '#e31a1c', '#fdbf6f', '#ff7f00', '#cab2d6', '#6a3d9a', '#ffff99', '#b15928']
# mpl.rcParams['axes.color_cycle'] = [list(clr) for clr in matplotlib.cm.brg(linspace(0,1,144))]

def bismark_analysis(input_cov_list, min_cov=10, fasta_fn=None, regions_fn=False, no_plot_scatters=False):
    """
    Run bismark comparison analysis - main function
    """
    # Load the data
    data = {}
    alldata = {}
    for fn in input_cov_list:
        alldata[fn], data[fn] = load_bismark_cov(fn, min_cov)
        if data is None:
            logging.error("Could not parse coverage file: {} Skipping..".format(fn))
            del (data[fn])
            del (alldata[fn])
            del (input_cov_list[fn])

    # Build nicer sample names
    sample_names = {}
    for fn in input_cov_list:
        name = os.path.basename(fn)
        remove = ['_val_1', '_val_2', '.fq', '.gz', '_bismark', '_bt2', '_bt1', '_pe', '_se', '.deduplicated', '.bismark', '.gwCov', '.cov']
        for r in remove:
            name = name.replace(r, '')
        name = name.rstrip('_1')
        name = name.rstrip('_2')
        sample_names[fn] = name

    # Plot dendrogram
    make_dendrogram(data, min_cov, sample_names)

    # Plot sample methylation histograms
    logging.info("Making histogram for combined samples.")
    plot_meth_histograms(data, sample_names)

    # Make one meth histogram per sample
    logging.info("Making one histogram per sample.")
    for fn in data.keys():
        plot_meth_histograms({fn: data[fn]}, {fn: sample_names[fn]}, output_fn="{}_histogram".format(sample_names[fn]))

    # Plot scatter plots and work out correlation scores
    if no_plot_scatters is not True:
        for i in range(0, len(input_cov_list)-1):
            for j in range(i+1, len(input_cov_list)):
                x = input_cov_list[j]
                y = input_cov_list[i]
                filenames = plot_meth_scatter(data[x], data[y], sample_names[x], sample_names[y], min_cov, no_plot_scatters)

    # TODO: Heatmap of sample correlation scores
    # TODO: Scatter plot of first two principal components

    # Do coverage analysis
    if fasta_fn is not None:
        if regions_fn is not None:
            regions = load_capture_regions(regions_fn)
        else:
            regions = None
        (fasta_ref, captured_cgs) = load_fasta_cpg(fasta_fn, regions)
        total_cg_count = len(fasta_ref)
        coverage_decay_plot(alldata, sample_names, total_cg_count, captured_cgs)


def load_bismark_cov(fn, min_cov=10):
    """
    Load a Bismark coverage file into memory. Takes normal coverage
    or genome-wide coverage files. The former are much, *much* faster
    to load.
    """
    logging.info("Loading {}".format(fn))
    alldata = defaultdict(dict)
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
                alldata[key]['coverage'] = cov
                alldata[key]['methylation'] = methp
                alldata[key]['chr'] = chrom
                alldata[key]['mb'] = int(float(c[1])/1000000)
                alldata[key]['pos'] = c[1]
                if cov >= min_cov:
                    data[key]['coverage'] = cov
                    data[key]['methylation'] = methp
                    data[key]['chr'] = chrom
                    data[key]['mb'] = int(float(c[1])/1000000)
                    data[key]['pos'] = c[1]

    except IOError as e:
        logging.error("Error loading coverage file: {}".format(fn))
        raise IOError(e)

    return alldata, data



def make_dendrogram(data, min_cov=10, sample_names=False, output_fn=False, plot_title=False):
    """
    Plots a hierarchical clustering dendrogram based on sample methylation
    scores (for positions seen in all samples)
    """

    # Output filename
    if output_fn is False:
        output_fn = "methylation_dendrogram"

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
    logging.info("Found {:.2e} ({:.0f}%) cytosines present in all datasets for dendrogram. Skipped {:.2e} ({:.0f}%) out of {:.2e} total observed (above coverage threshold of {}).".format(num_shared, percent_shared, skipped_keys, percent_skipped, num_total, min_cov))

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
    plt.figtext(0.99, 0.03, r'Calculated using {:.2e} cytosines found in all samples ({:.0f}% of all seen with coverage > {}X)'.format(num_shared, percent_shared, min_cov),
                horizontalalignment='right', fontsize=6, transform=axes.transAxes)

    # Save Plots
    png_fn = "{}.png".format(output_fn)
    pdf_fn = "{}.pdf".format(output_fn)
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)
    plt.close()

    # Return the filenames
    return {'png': png_fn, 'pdf': pdf_fn}


def plot_meth_histograms (data, sample_names, min_cov=10, output_fn=None):
    """
    Plots a histogram of methylation scores. Overlays all samples
    provided in a list
    """

    # Set up plot
    fig = plt.figure(figsize=(10,6))
    ax = plt.subplot(111)

    if len(data) == 1:
        onesample = True
    else:
        onesample = False

    # Go through datasets
    for fn, d in data.iteritems():

        # Get methylation scores
        meth = defaultdict(int)
        for p in d.values():
            m = "{:.0f}".format(float(p['methylation']))
            meth[m] += 1

        # Build plot list
        y_vals = []
        for m in range(0, 100):
            y_vals.append(meth[str(m)])

        # Plot
        if onesample is True:
            ax.bar(range(0, 100), y_vals, 1, color='#CCCCCC', linewidth=0.5)
            plt.title("Methylation Histogram: {} (Coverage > {}X)".format(sample_names[fn], min_cov))
        else:
            ax.plot(y_vals, linewidth=1)

    # Legend
    if onesample is False:
        # Give 10% room on the right hand side of the plot for the legend
        box = ax.get_position()
        ax.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        # Put a legend to the right of the current axis
        ax.legend(sample_names.values(), loc='upper left', bbox_to_anchor=(1, 1), fontsize=8)

    # Tidy plot
    plt.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)
    plt.xlabel('% Methylation')
    plt.ylabel('Number of CGs')
    if onesample is False:
        plt.title("Methylation Histogram (Coverage > {}X)".format(min_cov))

    # Save plot
    try:
        if output_fn is None:
            output_fn = "methylation_histogram"
        if onesample is False:
            plt.savefig("{}.png".format(output_fn))
            plt.savefig("{}.pdf".format(output_fn))
        else:
            if not os.path.exists("histograms/png"):
                os.makedirs("histograms/png")
            if not os.path.exists("histograms/pdf"):
                os.makedirs("histograms/pdf")
            plt.savefig("histograms/png/{}.png".format(output_fn))
            plt.savefig("histograms/pdf/{}.pdf".format(output_fn))
        plt.close()
    except IndexError:
        print(output_fn)



def plot_meth_scatter (sample_1, sample_2, x_lab, y_lab, min_cov=10, no_plot_scatters=False, output_fn=None):
    """
    Plots scatter plot of methylation score counts. Calculates r-squared
    correlation and spearman's correlation coefficient scores
    and prints these to STDOUT as well as including it in the graph.
    Returns a dict containing filenames of PNG and PDF graphs as a dict for
    each sample.
    """

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

    # Plot the scatter plot

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
    if not os.path.exists("scatter_plots/png"):
        os.makedirs("scatter_plots/png")
    if not os.path.exists("scatter_plots/pdf"):
        os.makedirs("scatter_plots/pdf")
    if output_fn is None:
        output_fn = "{}_{}_methylation_scatter".format(y_lab, x_lab)
    png_fn = "scatter_plots/png/{}.png".format(output_fn)
    pdf_fn = "scatter_plots/pdf/{}.pdf".format(output_fn)
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)
    plt.close()

    # Return the filenames and correlation scores
    return {'png': png_fn, 'pdf': pdf_fn}



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
    captured_cgs_count = 0
    if cap_regions is not None:
        captured_cgs = defaultdict(lambda:defaultdict(list))
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
            logging.info("  ..searching chromosome {} for CpGs - found {} so far, {} captured".format(chrom, len(fasta_ref), captured_cgs_count))

        # Look for +ve strand CpGs
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
                        mb = int(pos/1000000)
                        captured_cgs[chrom][mb].append(pos)
                        captured_cgs_count += 1

        # Look for -ve strand CpGs
        pos = sequence.find("GC")
        while pos >= 0:
            fasta_ref.append('{}_{}'.format(chrom, pos+1))
            pos = sequence.find("GC", pos + 1)
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
                        mb = int(pos/1000000)
                        captured_cgs[chrom][mb].append(pos)
                        captured_cgs_count += 1

    return (fasta_ref, captured_cgs)

@profile # Test with kernprof -l -v <command>
def coverage_decay_plot(data, sample_names=None, total_cg_count=False, captured_cgs=None):
    """
    Plots a coverage decay plot, y axis as percentage of genome-wide CpGs
    and y axis increasing coverage thresholds
    """
    # Set up variables for global plots
    g_vals = []
    if captured_cgs is not None:
        g_cap_vals = []

    # Collect counts of each coverage
    x_max = 200
    for fn, d in data.iteritems():

        # Labels
        try:
            sample_name = sample_names[fn]
        except KeyError:
            sample_name = fn

        logging.info("Creating coverage display plot for input: {}".format(sample_name))

        # Count the occurance of both types of coverage
        coverages = defaultdict(int)
        captured_coverages = defaultdict(int)
        for pos, arr in d.iteritems():
            if arr['coverage'] <= x_max:
                coverages[arr['coverage']] += 1
                if captured_cgs is not None:
                    if pos in captured_cgs[arr['chr']][arr['mb']]:
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
        g_vals.append([x_vals, y_vals])

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
            g_cap_vals.append([cap_x_vals, cap_y_vals])

        # Draw the plot
        fig = plt.figure(figsize=(10,6))
        # Bar plot if only one dataset
        if captured_cgs is None:
            ax = plt.subplot(111)
            ax.bar(x_vals, y_vals, 1)
        # Lines if plotting two things
        else:
            plt.plot(x_vals, y_vals, 'b-')
            plt.plot(cap_x_vals, cap_y_vals, 'g--')

        # Tidy axes
        plt.xlim(0,x_max)
        plt.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)

        plt.xlabel('CpG Coverage')
        plt.ylabel('Percentage of CGs')
        plt.title("Coverage Decay Plot: {}".format(sample_name))

        # Save the individual plot
        if not os.path.exists("coverage/png"):
            os.makedirs("coverage/png")
        if not os.path.exists("coverage/pdf"):
            os.makedirs("coverage/pdf")
        output_fn = "{}_coverage_decay".format(sample_name)
        plt.savefig('coverage/png/{}.png'.format(output_fn))
        plt.savefig('coverage/pdf/{}.pdf'.format(output_fn))
        plt.close()

    # Draw the global plot
    fig = plt.figure(figsize=(10,6))
    for vals in g_vals:
        plt.plot(vals[0], vals[1])

    # Tidy axes
    plt.xlim(0,x_max)
    plt.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)

    plt.xlabel('CpG Coverage')
    plt.ylabel('Percentage of CGs')
    plt.title("Coverage Decay Plot")

    # Save the individual plot
    output_fn = "all_coverage_decay"
    plt.savefig('{}.png'.format(output_fn))
    plt.savefig('{}.pdf'.format(output_fn))
    plt.close()

    # Draw the global captured plot
    if captured_cgs is not None:
        fig = plt.figure(figsize=(10,6))
        for vals in g_cap_vals:
            plt.plot(vals[0], vals[1])

        # Tidy axes
        plt.xlim(0,x_max)
        plt.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)

        plt.xlabel('CpG Coverage')
        plt.ylabel('Percentage of CGs')
        plt.title("Coverage Decay Plot - Captured Regions")

        # Save the individual plot
        output_fn = "capture_coverage_decay"
        plt.savefig('{}.png'.format(output_fn))
        plt.savefig('{}.pdf'.format(output_fn))
        plt.close()


if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("Compare bismark methylation coverage data between samples")
    parser.add_argument("-m", "--min-cov", dest="min_cov", default=10,
                        help="Minimum coverage to use for scatter plots and dendrogram")
    parser.add_argument("-f", "--fasta", dest="fasta_fn",
                        help="Genome Fasta reference for coverage analysis")
    parser.add_argument("-r", "--regions", dest="regions_fn",
                        help="BED file containing capture regions")
    parser.add_argument("-s", "--no-scatter", dest="no_plot_scatters", action='store_true',
                        help="Don't plot pairwise scatter plots.")
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
