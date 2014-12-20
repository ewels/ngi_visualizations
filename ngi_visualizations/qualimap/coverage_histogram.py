#!/usr/bin/python
"""
coverage_histogram.py

Takes results from the Qualimap coverage_histogram.txt results file
and plots a nice looking graph using matplotlib
"""

from __future__ import print_function

import argparse
import logging
import os

# Import matplot lib but avoid default X environment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def plot_coverage_histogram (coverage_data, output_fn='coverage', min_x='0', max_x='Auto'):
    """
    Main function. Takes input file and makes a plot.
    """
    # Sort out the incoming variables
    try:
        min_x = int(min_x)
    except (ValueError, TypeError):
        min_x = None
    try:
        max_x = int(max_x)
    except (ValueError, TypeError):
        if max_x == 'None':
            max_x = None
        else:
            max_x = False
    
    # Load in the data
    fn = os.path.realpath(coverage_data)
    counts={}
    try:
        with open(fn, 'r') as fh:
            next(fh) # skip the header
            for line in fh:
                (coverage, count) = line.split(None, 1)
                coverage = int(round(float(coverage)))
                count = float(count)
                counts[coverage] = count
                

    except IOError as e:
        logging.error("Could not load input file: {}".format(fn))
        raise
    
    # Find mean
    num_counts = sum(counts.values())
    mean_coverage = sum([cou * cov for (cou, cov) in counts.iteritems()]) / num_counts
    
    # Find median
    cum_counts = 0
    median_coverage = None
    for thiscov, thiscount in counts.iteritems():
        cum_counts += thiscount
        if cum_counts >= num_counts/2:
            median_coverage = thiscov
            break
    
    # Set max_x automatically if we want to
    if max_x == False and median_coverage is not None:
        max_x = median_coverage * 2
    
    # Collect the x and y points that we need and ignore the ones that we don't
    x = []
    y = []
    for thiscov, thiscount in counts.iteritems():
        if (min_x is None or thiscov >= min_x) and (max_x is None or thiscov <= max_x):
            x.append(thiscov)
            y.append(thiscount)
        elif max_x is not None and thiscov > max_x:
            break
    
    # Check that we found something
    if len(y) == 0:
        raise EOFError ("Unable to find any data in input file within coverage range")
    
    # Set up the plot
    fig = plt.figure(figsize=(8,3.4), tight_layout={'rect':(0,0.04,1,1)})
    axes = fig.add_subplot(111)
    [i.set_linewidth(0.5) for i in axes.spines.itervalues()] # thinner border
    
    # Plot
    axes.bar(x, y, width=1, linewidth=0.5, color="#b3cde3")
    
    # Plot the mean as a dashed line
    axes.axvline(mean_coverage, color='black', linestyle='--', linewidth=1)

    # Set axis limit if we need to
    if min_x is not None and max_x is not None:
        axes.set_xlim([min_x, max_x])
    elif min_x is not None:
        axes.set_xlim(left=int(min_x))
    elif max_x is not None:
        axes.set_xlim(right=max_x)

    # Tidy the axes
    axes.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)
    axes.grid(True, zorder=0, which='both', axis='y', linestyle='-', color='#EDEDED', linewidth=1)
    axes.set_axisbelow(True)
    
    # Make the x axis labels have an X
    axes.set_xticklabels(["%dX" % d for d in axes.get_xticks()])
    
    # Labels
    matplotlib.rcParams['mathtext.default'] = 'regular'
    plt.xlabel(r"Coverage")
    plt.ylabel(r'Genome Bin Counts')
    plt.text(0.5, -0.25, 'Mean Coverage: {0:.2f}X'.format(mean_coverage),
                horizontalalignment='center', fontsize=8, transform = axes.transAxes)

    # SAVE OUTPUT
    png_fn = "{}.png".format(output_fn)
    pdf_fn = "{}.pdf".format(output_fn)
    logging.info("Saving to {} and {}".format(png_fn, pdf_fn))
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)
    
    # Close the plot
    plt.close(fig)



if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("coverage_histogram.py", description="Plot coverage histogram")
    parser.add_argument("-o", "--output", dest="output_fn", default='coverage',
                        help="Plot output filename base. Default: coverage.png / .pdf")
    parser.add_argument("-x", "--min_x", dest="min_x", default='0',
                        help="Minimum x axis limit. Use 'None' for data limit. Default: 0")
    parser.add_argument("-m", "--max_x", dest="max_x", default='Auto',
                            help="Maximum x axis limit. Use 'None' for data limit or 'Auto' for 2 * median. Default: Auto")
    parser.add_argument("-l", "--log", dest="log_level", default='info', choices=['debug', 'info', 'warning'],
                        help="Level of log messages to display")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("coverage_data", metavar='<coverage data file>',
                        help="Data input - usually raw_data_qualimapReport/coverage_histogram.txt")
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
    
    # Call plot_observed_genes()
    plot_coverage_histogram(**kwargs)