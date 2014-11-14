#!/usr/bin/python
"""
snpEff_plots.py

Takes results from the snpEff snpEff_summary.csv results file
and plots some nice looking graphs using matplotlib
"""

from __future__ import print_function

import argparse
import logging
import numpy
import os

# Import matplot lib but avoid default X environment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_snpEff (summary_fn, output_fn='effect_types', logx=False):
    """
    Main function. Takes summary input file and makes some plots.
    """    
    # Load in the data
    fn = os.path.realpath(summary_fn)
    types = {}
    counts = {}
    percents = {}
    try:
        with open(fn, 'r') as fh:
            recording = False
            for line in fh:
                line = line.strip()
                if recording is False and line == '# Count by effects':
                    recording = True
                elif recording and line[:1] == '#':
                    recording = False
                if recording:
                    if line.count(',') >= 2:
                        (effect_type, count, percent) = line.split(' , ', 2)
                        if count.isdigit():
                            effect_type_nice = effect_type.replace('_', ' ').title().replace('Utr', 'UTR')
                            count = float(count)
                            percent = float(percent[:-1])
                            types[effect_type] = effect_type_nice
                            counts[effect_type] = count
                            percents[effect_type] = percent
    except IOError as e:
        logging.error("Could not load input file: {}".format(fn))
        raise IOError(e)
    
    # Check that we found something
    if len(counts) == 0:
        raise EOFError ("Unable to find any effect type counts in input file")
    
    # Prepare the sorted values
    plt_labels = []
    plt_counts = []
    percent_factor = 0
    for effect in sorted(counts, key=counts.get):
        plt_labels.append(types[effect])
        plt_counts.append(counts[effect])
        percent_factor = percents[effect] / counts[effect]
    
    # Set up the plot
    fig = plt.figure(figsize=(8,5), tight_layout={'rect':(0,0.04,1,0.95)})
    axes = fig.add_subplot(111)
    [i.set_linewidth(0.5) for i in axes.spines.itervalues()] # thinner border
    
    # Plot
    min_x = 0
    if logx:
        min_x = 1
    ypos = numpy.arange(1, len(plt_labels)+1)
    axes.bar(min_x, 0.8, plt_counts, ypos, log=logx, orientation='horizontal', linewidth=0.5, color="#fbb4ae")

    # Y axis
    axes.set_yticks(ypos+0.3)
    axes.set_yticklabels(plt_labels)
    axes.tick_params(left=False, right=False)
    
    # X axis
    if logx:
        axes.set_xscale('log')
    axes.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)
    axes.grid(True, zorder=0, which='both', axis='x', linestyle='-', color='#EDEDED', linewidth=1)
    axes.set_axisbelow(True)
    
    # Make the secondary percentage axis
    ax2 = axes.twiny()
    ax2.tick_params(which='both', labelsize=8, direction='out', bottom=False, right=False)
    if logx:
        ax2.set_xscale('log')
    ax2.set_xlim(axes.get_xlim())
    ax1_ticks = axes.get_xticks()
    # I have no idea why I have to get rid of these two elements....
    ax1_ticks = ax1_ticks[1:-1]
    ax2.set_xticks(ax1_ticks)
    
    # Set the secondary axis percentage labels
    def percent_total(counts):
        y = [x*percent_factor for x in counts]
        return ["%.2f%%" % z for z in y]
    ax2_labels = percent_total(ax2.get_xticks())
    ax2.set_xticklabels(ax2_labels)
    
    # Labels
    plt.text(0.5, 1.1, 'SNP Effects', horizontalalignment='center',
                    fontsize=16, transform=axes.transAxes)
    matplotlib.rcParams['mathtext.default'] = 'regular'
    axes.set_xlabel(r"Number of Effects")
    plt.ylabel(r'Effect Type')

    # SAVE OUTPUT
    png_fn = "{}.png".format(output_fn)
    pdf_fn = "{}.pdf".format(output_fn)
    logging.info("Saving to {} and {}".format(png_fn, pdf_fn))
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)




if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("snpEff_plots.py", description="Plot snpEFF graphs")
    parser.add_argument("-o", "--output", dest="output_fn", default='effect_types',
                        help="Plot output filename base. Default: effect_types.png / .pdf")
    parser.add_argument("-x", "--logx", dest="logx", type=bool, default=False,
                        help="Use a log scale on the x axis?")
    parser.add_argument("-l", "--log", dest="log_level", default='info', choices=['debug', 'info', 'warning'],
                        help="Level of log messages to display")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("summary_fn", metavar='<snpEff summary file>',
                        help="Data input - usually snpEff_summary.csv")
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
    plot_snpEff(**kwargs)