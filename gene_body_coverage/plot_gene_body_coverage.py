#!/usr/bin/python
"""
plot_gene_body_coverage.py

Takes results from the RSeQC geneBody_coverage.py script and plots them
using matplotlib
"""

from __future__ import print_function

import argparse
import logging
import os

# Import matplot lib but avoid default X environment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_gene_body_coverage (coverage_files, output_fn='geneBodyCoverage'):
    """
    Main function. Takes input files and makes a plot.
    """
    
    # Set up the plot
    fig = plt.figure()
    axes = fig.add_subplot(111)
    plt.subplots_adjust(right=0.7)
    colours = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']*10
    i = 0
    for fn in coverage_files:
        fn = os.path.realpath(fn)
        values = [0] * 101
        try:
            with open(fn, 'r') as fh:
                for line in fh:
                    (percentile, count) = line.split(None, 1)
                    if percentile.isdigit() is False:
                        continue
                    values[int(percentile)] = float(count) / 1000000
        except IOError as e:
            logging.error("Could not load input file: {}".format(fn))
            raise IOError(e)
        
        sample_name = os.path.splitext(os.path.splitext(os.path.basename(fn))[0])[0]
        axes.plot(values, label=sample_name, color=colours[i])
        i += 1
    
    # Tidy the axes
    axes.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)
    
    # Make the x axis labels percentages
    axes.set_xticklabels(["%d%%" % d for d in axes.get_xticks()])
    
    # Labels
    matplotlib.rcParams['mathtext.default'] = 'regular'
    plt.xlabel(r"Gene body position ($5' \rightarrow 3'$)")
    plt.ylabel(r'Read Count ($\times 10^6$)')
    plt.title('Gene Body Coverage')
    
    # Legend
    axes.legend(loc='upper left', bbox_to_anchor = (1.02, 1.02), fontsize=8)
    
    # SAVE OUTPUT
    png_fn = "{}.png".format(output_fn)
    pdf_fn = "{}.pdf".format(output_fn)
    logging.info("Saving to {} and {}".format(png_fn, pdf_fn))
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)
    
    
if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("plot_gene_body_coverage.py", description="Plot number of observed genes at increasing read depths")
    parser.add_argument("-o", "--output", dest="output_fn", default='geneBodyCoverage',
                        help="Plot output filename base. Default: geneBodyCoverage.png / .pdf")
    parser.add_argument("-l", "--log", dest="log_level", default='info', choices=['debug', 'info', 'warning'],
                        help="Level of log messages to display")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("coverage_files", metavar='<gene body coverage files>', nargs="+",
                        help="List of gene body coverage results files")
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
    plot_gene_body_coverage(**kwargs)