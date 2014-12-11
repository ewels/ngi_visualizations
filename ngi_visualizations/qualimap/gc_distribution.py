#!/usr/bin/python
"""
gc_distribution.py

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


def plot_genome_fraction_coverage (gc_distribution_input, output_fn='gc_distribution', min_x='None', max_x='None', picard=None, reference_label='Reference Genome (hg19)'):
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
        max_x = None
    
    # Load in the data
    fn = os.path.realpath(gc_distribution_input)
    x1 = [0]
    y1 = [0]
    x2 = [0]
    y2 = [0]
    try:
        with open(fn, 'r') as fh:
            next(fh) # skip the header
            for line in fh:
                sections = line.split(None, 2)
                gc = int(round(float(sections[0])))
                percentage = float(sections[1]) * 100
                try:
                    reference = float(sections[2]) * 100
                except IndexError:
                    reference = False
                if (min_x is None or gc >= min_x) and (max_x is None or gc <= max_x):
                    x1.append(gc)
                    y1.append(percentage)
                    if reference:
                        x2.append(gc)
                        y2.append(reference)

    except IOError as e:
        logging.error("Could not load input file: {}".format(fn))
        raise IOError(e)
    
    # Check that we found something
    if len(y1) == 1:
        raise EOFError ("Unable to find any data in input file")
    
    # Load picard data if we have it
    x3 = [0]
    y3 = [0]
    if picard is not None:
        fn = os.path.realpath(picard)
        try:
            with open(fn, 'r') as fh:
                for line in fh:
                    line = line.strip()
                    if len(line) == 0 or line[:1] == '#':
                        continue
                    sections = line.split(None, 6)
                    if sections[0][:-1].isdigit() is False:
                        continue
                    gc = int(round(float(sections[0])))
                    coverage = float(sections[4])
                    if (min_x is None or gc >= min_x) and (max_x is None or gc <= max_x):
                        x3.append(gc)
                        y3.append(coverage)
        except IOError as e:
            logging.error("Could not load picard CollectGcBiasMetrics file: {}".format(fn))
            logging.error("Skipping normalised coverage plotting")
            pass
    
    # Set up the plot
    fig = plt.figure(figsize=(8,3.4), tight_layout={'rect':(0,0,1,1)})
    axes = fig.add_subplot(111)
    [i.set_linewidth(0.5) for i in axes.spines.itervalues()] # thinner border
    
    # Plot the histogram
    axes.bar(x1, y1, width=1, linewidth=0.5, color="#decbe4")
    
    # Plot the reference, if we have it
    lns = []
    if len(y2) > 1:
        lns += axes.plot(x2, y2, 'r--', label=reference_label)
    
    # Plot the normalised coverage, if we have it
    if len(y3) > 1:
        ax2 = axes.twinx()
        lns += ax2.plot(x3, y3, 'b--', label="Normalised Coverage")
        ax2.set_ylim([0,2])
        ax2.axhline(1, color='black', linestyle=':', linewidth=1)
    
    # Plot a legend if we need it
    if len(y2) > 1 or len(y3) > 1:
        labs = [l.get_label() for l in lns]
        if len(y3) > 1:
            leg = ax2.legend(lns, labs, prop={'size':8}, handlelength=2.5)
        else:
            leg = axes.legend(lns, labs, prop={'size':8}, handlelength=2.5)
        leg.get_frame().set_linewidth(0.0)

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
    
    # Now for second y axis if we have it
    if len(y3) > 1:
        ax2.tick_params(which='both', labelsize=8, direction='out', top=False, left=False)

    # Make the x and y axis labels have a %
    axes.set_yticklabels(["%.1f%%" % d for d in axes.get_yticks()])
    axes.set_xticklabels(["%d%%" % d for d in axes.get_xticks()])
    
    # Labels
    matplotlib.rcParams['mathtext.default'] = 'regular'
    axes.set_xlabel(r"% GC Content")
    axes.set_ylabel(r'Percentage of Library')
    if len(y3) > 1:
        ax2.set_ylabel(r'Normalised Coverage')
    
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
    parser = argparse.ArgumentParser("gc_distribution.py", description="Plot GC-content Distribution")
    parser.add_argument("-o", "--output", dest="output_fn", default='gc_distribution',
                        help="Plot output filename base. Default: gc_distribution.png / .pdf")
    parser.add_argument("-x", "--min_x", dest="min_x", default='0',
                        help="Minimum x axis limit. Use 'None' for data limit. Default: None")
    parser.add_argument("-m", "--max_x", dest="max_x", default='100',
                            help="Maximum x axis limit. Use 'None' for data limit. Default: None")
    parser.add_argument("-p", "--picard", dest="picard", default=None,
                            help="Picard CollectGcBiasMetrics output file to plot normalised GC content. Default: None")                            
    parser.add_argument("-r", "--ref_label", dest="reference_label", default='Reference Genome (hg19)',
                            help="Legend for the reference data if present. Default: 'Reference Genome (hg19)'")                        
    parser.add_argument("-l", "--log", dest="log_level", default='info', choices=['debug', 'info', 'warning'],
                        help="Level of log messages to display")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("gc_distribution_input", metavar='<gc content distribution data file>',
                        help="Data input - usually raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt")
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