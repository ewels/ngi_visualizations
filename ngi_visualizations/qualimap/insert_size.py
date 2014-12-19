#!/usr/bin/python
"""
insert_size.py

Takes results from the Qualimap insert_size_histogram.txt results file
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


def plot_insert_size_histogram (insertsize_data, output_fn='insert_size', min_x='None', max_x='Auto', bin_size=10):
    """
    Main function. Takes input file and makes a plot.
    """
    # Sort out the incoming variables
    try:
        min_x = int(min_x)
    except ValueError, TypeError:
        min_x = None
    try:
        max_x = int(max_x)
    except ValueError, TypeError:
        if max_x == 'None':
            max_x = None
        else:
            max_x = False
        
    # Load in the data
    fn = os.path.realpath(insertsize_data)
    counts = {}
    try:
        with open(fn, 'r') as fh:
            next(fh) # skip the header
            for line in fh:
                (insertsize, count) = line.split(None, 1)
                insertsize = int(round(float(insertsize)))
                count = float(count) / 1000000
                counts[insertsize] = count

    except IOError as e:
        logging.error("Could not load input file: {}".format(fn))
        raise IOError(e)
    
    # Find mean
    num_counts = sum(counts.values())
    mean_insert_size = sum([ins * cov for (ins, cov) in counts.iteritems()]) / num_counts
    
    # Find median
    cum_counts = 0
    median_insert_size = None
    for thisins, thiscount in counts.iteritems():
        cum_counts += thiscount
        if cum_counts >= num_counts/2:
            median_insert_size = thisins
            break
    
    # Set max_x automatically if we want to
    if max_x == False and median_insert_size is not None:
        max_x = median_insert_size * 2
    
    # Collect the x and y points that we need and ignore the ones that we don't
    x = []
    y = []
    bincounts = 0
    for thisins, thiscount in counts.iteritems():
        if (min_x is None or thisins >= min_x) and (max_x is None or thisins <= max_x):
            bincounts += thiscount
            if thisins % bin_size == 0:
                x.append(thisins)
                y.append(bincounts)
                bincounts = 0
        elif max_x is not None and thisins > max_x:
            break
    
    # Check that we found something
    if len(y) == 0:
        raise EOFError ("Unable to find any data in input file")
    
    # Set up the plot
    fig = plt.figure(figsize=(8,3.4), tight_layout={'rect':(0,0.04,1,1)})
    axes = fig.add_subplot(111)
    [i.set_linewidth(0.5) for i in axes.spines.itervalues()] # thinner border
    
    # Plot
    axes.bar(x, y, width=bin_size, linewidth=0.5, color="#ccebc5")
    
    # Plot the mean as a dashed line
    axes.axvline(mean_insert_size, color='black', linestyle='--', linewidth=1)

    # Set axis limit if we need to
    if min_x and max_x:
        axes.set_xlim([min_x, max_x])
    elif min_x:
        axes.set_xlim(left=min_x)
    elif max_x:
        axes.set_xlim(right=max_x)

    # Tidy the axes
    axes.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)
    axes.grid(True, zorder=0, which='both', axis='y', linestyle='-', color='#EDEDED', linewidth=1)
    axes.set_axisbelow(True)
    # axes.set_xlim([0,60])
    
    # Make the x axis labels have an X
    axes.set_xticklabels(["%d bp" % d for d in axes.get_xticks()])
    
    # Labels
    matplotlib.rcParams['mathtext.default'] = 'regular'
    plt.xlabel(r"Insert Size (bp)")
    plt.ylabel(r'Number of Reads ($\times 10^6$)')
    plt.text(0.5, -0.25, 'Mean Insert Size: {:.0f} bp'.format(mean_insert_size),
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
    parser = argparse.ArgumentParser("insert_size.py", description="Plot insertsize histogram")
    parser.add_argument("-o", "--output", dest="output_fn", default='insertsize',
                        help="Plot output filename base. Default: insert_size.png / .pdf")
    parser.add_argument("-x", "--min_x", dest="min_x", default='None',
                        help="Minimum x axis limit. Use 'None' for data limit. Default: None")
    parser.add_argument("-m", "--max_x", dest="max_x", default='Auto',
                            help="Maximum x axis limit. Use 'None' for data limit or 'Auto' for 2 * median. Default: Auto")
    parser.add_argument("-b", "--bin_size", dest="bin_size", type=int, default=10,
                            help="Histogram bin size. Default: 10 bp")
    parser.add_argument("-l", "--log", dest="log_level", default='info', choices=['debug', 'info', 'warning'],
                        help="Level of log messages to display")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("insertsize_data", metavar='<insert size data file>',
                        help="Data input - usually raw_data_qualimapReport/insert_size_histogram.txt")
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
    plot_insert_size_histogram(**kwargs)