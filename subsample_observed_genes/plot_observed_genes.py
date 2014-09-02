#!/usr/bin/python
"""
plot_observed_genes.py

Takes results from cufflinks for a range of subsampled BAM files
and plots the number of genes with FPKM > 1 at increasing sampling
rates.

Input must be a space separated list of directories. They must be
named with the structure:
  */<original_bam_filename>_<subsample_fraction>/
eg:
  path/to/data/inputFile1.bam_0.1/
  path/to/data/inputFile1.bam_0.2/
  ..
  path/to/data/inputFile1.bam_1.0/
"""

from __future__ import print_function

import argparse
from collections import defaultdict
import logging
import numpy
import os
import re

# Import matplot lib but avoid default X environment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def plot_observed_genes (input_dirs, output_fn='gene_counts.png', counts_fn=None):
    """
    Main function. Takes input files and makes a plot.
    """
    # Parse the directory names
    samples = parse_input_dirnames(input_dirs)
    
    # Find the counts
    counts = count_cufflinks_observed_genes(samples)
    
    # Plot the counts
    filenames = plot_gene_counts(counts)
    
    
    
def parse_input_dirnames (input_dirs):
    """
    Take a list of input directories and sort them into sample and proportion.
    Returns a dict with keys = sample names. Each var is a dict with keys
    equal to the proportion number and value the path to the genes.fpkm_tracking
    If we find a parsing error then we die.
    """
    # Check that we have enough directories
    if len(input_dirs) < 2:
        raise IOError("Fatal error - we need at least two input directories!")
    
    # Go through directories
    samples = defaultdict(lambda: defaultdict(str))
    re_pattern = re.compile('^(.*)_([\d\.]+)$')
    for dirname in input_dirs:
        # check that we're a dir and that we exists
        if not os.path.isdir(dirname):
            if os.path.exists(dirname):
                
                raise IOError("Fatal error - not a directory: {}".format(dirname))
            else:
                raise IOError("Fatal error - can't find input directory: {}".format(dirname))
        # Strip the basename from the path
        dirbase = os.path.basename(dirname)
        # Parse the names and assign to the array
        m = re_pattern.match(dirbase)
        if m is None:
            raise Exception("Fatal error - couldn't match the directory name {}".format(dirbase))
        else:
            try:
                sample = m.group(1)
                proportion = m.group(2)
                filename = '{}/genes.fpkm_tracking'.format(os.path.realpath(dirname))
                samples[sample][proportion] = filename
            except IndexError as e:
                logging.error("Fatal error - problem with the directory name matches for {} - {}".
                            format(dirname, e))
                raise IndexError(e)
    
    # Check that we have some samples
    if len(samples) == 0:
        raise Exception("Error - couldn't find any input samples!")
    
    # Return our parsed dict!
    logging.info("Found {} samples..".format(len(samples)))
    return samples
    
def count_cufflinks_observed_genes (samples):
    """
    Reads cufflinks genes.fpkm_tracking files and counts the number of
    genes with FPKM > 1
    
    Input: Nested dict of structure samples[sample_fn][proportion] = fn
    where fn is the path to genes.fpkm_tracking
    
    Output: Nested dict with same structure, but values as gene count
    """
    
    counts = defaultdict(lambda: defaultdict(int))
    for sample in samples:
        logging.info("Processing {}".format(sample))
        for proportion in sorted(samples[sample]):
            logging.debug("  ..proportion {}".format(proportion))
            counts[sample][proportion] = 0
            file = samples[sample][proportion]
            if os.path.isfile(file):
                try:
                    with open(file, 'r') as fh:
                        # Skip the header
                        next(fh)
                        # Iterate through the file
                        for line in fh:
                            line = line.strip()
                            cols = line.split("\t")
                            FPKM = cols[9]
                            if float(FPKM) > 1:
                                counts[sample][proportion] += 1
                except IOError as e:
                    logging.error("Error loading cufflinks input file: {}".format(file))
                    raise IOError(e)
            else:
                raise IOError("Couldn't find cufflinks input file: {}".format(file))
            
            logging.info("{} - {} = {}".format(sample, proportion, counts[sample][proportion]))
    
    return counts



def plot_gene_counts (counts, title="Subsampled Gene Discovery Rates"):
    """
    Plots line graph of numbers of genes with FPKM > 1 at increasing subsampling
    levels, using matplotlib pyplot
    Input: counts[sample_id][subsample_proportion] = int <gene_counts>
    Input: title=Plot title
    Returns a dict containing filenames of PNG and PDF graphs as a dict for 
    each sample. Filenames are generated from sample IDs
    """
    
    # SET UP PLOT
    fig = plt.figure()
    axes = fig.add_subplot(111)
    plt.subplots_adjust(right=0.7)
    colours = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']*10
    
    i = 0
    for sample in counts:
        
        proportions = [0]
        gene_counts = [0]
        
        for p in sorted(counts[sample]):
            gene_counts.append(counts[sample][p])
            proportions.append(p)
        
        axes.plot(proportions, gene_counts, color=colours[i], label=sample)
        i += 1
    
    # Tidy axes
    axes.set_xlim(0,1)
    axes.grid(True, zorder=0, which='both', axis='y', linestyle='-', color='#EDEDED', linewidth=1)
    axes.set_axisbelow(True)
    axes.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)
    
    # Labels
    plt.xlabel('Proportion of sample')
    plt.ylabel('Number of genes with FPKM > 1')
    plt.title(title)

    # Legend
    axes.legend(loc='upper left', bbox_to_anchor = (1.02, 1.02), fontsize=8)
    
    # SAVE OUTPUT
    png_fn = "GeneCounts.png".format(sample)
    pdf_fn = "GeneCounts.pdf".format(sample)
    logging.info("Saving to {} and {}".format(png_fn, pdf_fn))
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)

    # Return the filenames
    return {'png': png_fn, 'pdf': pdf_fn}


    
    
if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("Plot number of observed genes at increasing read depths")
    parser.add_argument("-c", "--counts", dest="counts_fn", default=None,
                        help="Filename of counts file to use for x axis instead of proportions. See README for further info.")
    parser.add_argument("-o", "--output", dest="output_fn", default='gene_counts.png',
                        help="Plot output filename. Default: gene_counts.png")
    parser.add_argument("-l", "--log", dest="log_level", default='info', choices=['debug', 'info', 'warning'],
                        help="Level of log messages to display")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("input_dirs", metavar='<cufflinks directory>', nargs="+",
                        help="List of cufflinks results directories")
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
    plot_observed_genes(**kwargs)