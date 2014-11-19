#!/usr/bin/python
"""
genome_fraction_coverage.py

Takes results from the Qualimap genome_fraction_coverage.txt results file
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


def plot_genome_fraction_coverage (fraction_data, output_fn='genome_fraction', min_x='None', max_x='None'):
    """
    Main function. Takes input file and makes a plot.
    """
    # Sort out the incoming variables
    if min_x.isdigit(): min_x = int(min_x)
    else: min_x = None
    if max_x.isdigit(): max_x = int(max_x)
    else: max_x = None
    
    # Load in the data
    fn = os.path.realpath(fraction_data)
    x = [0]
    y = [100]
    eighty_pc_coverage = 0
    thirty_x_pc = 100
    try:
        with open(fn, 'r') as fh:
            next(fh) # skip the header
            for line in fh:
                (coverage, percentage) = line.split(None, 1)
                coverage = int(round(float(coverage)))
                percentage = float(percentage)
                if (min_x is None or coverage >= min_x) and (max_x is None or coverage <= max_x):
                    x.append(coverage)
                    y.append(percentage)
                if percentage >= 80 and eighty_pc_coverage < coverage:
                    eighty_pc_coverage = coverage
                if coverage <= 30 and thirty_x_pc > percentage:
                    thirty_x_pc = percentage

    except IOError as e:
        logging.error("Could not load input file: {}".format(fn))
        raise IOError(e)
    
    # Check that we found something
    if len(y) == 1:
        raise EOFError ("Unable to find any data in input file")
    
    # Set up the plot
    fig = plt.figure(figsize=(8,3.4), tight_layout={'rect':(0,0.04,1,1)})
    axes = fig.add_subplot(111)
    [i.set_linewidth(0.5) for i in axes.spines.itervalues()] # thinner border
    
    # Plot
    axes.bar(x, y, width=1, linewidth=0.5, color="#fbb4ae")

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
    
    # Make the y axis labels have a %
    axes.set_yticklabels(["%d%%" % d for d in axes.get_yticks()])
    
    # Make the x axis labels have an X
    axes.set_xticklabels(["%dX" % d for d in axes.get_xticks()])
    
    # Labels
    matplotlib.rcParams['mathtext.default'] = 'regular'
    plt.xlabel(r"Coverage")
    plt.ylabel(r'% Ref. Genome $\geq$ Coverage')
    plt.text(0.5, -0.25, r'80% of reference has $\geq$ {:.0f}X coverage. {:.2f}% of reference has $\geq$ 30X coverage.'.format(eighty_pc_coverage, thirty_x_pc),
                horizontalalignment='center', fontsize=8, transform = axes.transAxes)
    
    # SAVE OUTPUT
    png_fn = "{}.png".format(output_fn)
    pdf_fn = "{}.pdf".format(output_fn)
    logging.info("Saving to {} and {}".format(png_fn, pdf_fn))
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)




if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("genome_fraction_coverage.py", description="Plot genome fraction coverage")
    parser.add_argument("-o", "--output", dest="output_fn", default='genome_fraction',
                        help="Plot output filename base. Default: genome_fraction.png / .pdf")
    parser.add_argument("-x", "--min_x", dest="min_x", default='None',
                        help="Minimum x axis limit. Use 'None' for data limit. Default: None")
    parser.add_argument("-m", "--max_x", dest="max_x", default='None',
                            help="Maximum x axis limit. Use 'None' for data limit. Default: None")
    parser.add_argument("-l", "--log", dest="log_level", default='info', choices=['debug', 'info', 'warning'],
                        help="Level of log messages to display")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("fraction_data", metavar='<genome fraction data file>',
                        help="Data input - usually raw_data_qualimapReport/genome_fraction_coverage.txt")
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
    plot_genome_fraction_coverage(**kwargs)