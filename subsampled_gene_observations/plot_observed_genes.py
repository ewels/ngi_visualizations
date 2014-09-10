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

def plot_observed_genes (input_dirs, fpkm_cutoff=0, output_fn='gene_counts', read_counts_fn=None):
    """
    Main function. Takes input files and makes a plot.
    """
    # Parse the directory names
    samples = parse_input_dirnames(input_dirs)
    
    # Find the gene counts
    gene_counts = count_cufflinks_observed_genes(samples, fpkm_cutoff)
    
    # Try to find read counts if we have the cufflinks logs
    if read_counts_fn is not None:
        read_counts = get_read_counts(read_counts_fn, samples)
    else:
        read_counts = None
    
    # Plot the proportions
    filenames = plot_gene_counts(gene_counts, fpkm_cutoff, None, "{}_proportions".format(output_fn))
    
    # Plot the counts if we have them
    if read_counts is not None:
        proportion_filenames = plot_gene_counts(gene_counts, fpkm_cutoff, read_counts, output_fn)
    

def parse_filename (fn):
    """
    Helper function - takes a file name as input and returns a dict with 
    the sample name and proportion.
    """
    re_pattern = re.compile('^(.*)_([\d\.]+)$')
    m = re_pattern.match(fn)
    if m is None:
        return None
    else:
        sample = m.group(1)
        proportion = m.group(2)
        match = {"sample":sample, "proportion":proportion}
        return match
    
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
    for dirname in sorted(input_dirs):
        # check that we're a dir and that we exists
        if not os.path.isdir(dirname):
            raise IOError("Fatal error - can't find input directory: {}".format(dirname))
        # Try the pattern matching
        try:
            filename = '{}/genes.fpkm_tracking'.format(os.path.realpath(dirname))
            match = parse_filename(os.path.basename(dirname))
            samples[match['sample']][match['proportion']] = filename
        except IndexError:
            raise IOError("Fatal error - couldn't match the directory name {}".format(dirbase))
    
    # Check that we have some samples
    if len(samples) == 0:
        raise Exception("Error - couldn't find any input samples!")
    
    # Return our parsed dict!
    logging.info("Found {} samples..".format(len(samples)))
    return samples
    
def count_cufflinks_observed_genes (samples, fpkm_cutoff=0):
    """
    Reads cufflinks genes.fpkm_tracking files and counts the number of
    genes with FPKM > 1
    
    Input: Nested dict of structure samples[sample_fn][proportion] = fn
    where fn is the path to genes.fpkm_tracking
    
    Output: Nested dict with same structure, but values as gene count
    """
    
    fpkm_cutoff = float(fpkm_cutoff)
    
    counts = defaultdict(lambda: defaultdict(int))
    for sample in sorted(samples):
        logging.info("Processing {}".format(sample))
        for proportion in sorted(samples[sample]):
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
                            if float(FPKM) > fpkm_cutoff:
                                counts[sample][proportion] += 1
                except IOError as e:
                    logging.error("Error loading cufflinks input file: {}".format(file))
                    raise IOError(e)
            else:
                raise IOError("Couldn't find cufflinks input file: {}".format(file))
            
            logging.info("  ..{} = {}".format(proportion, counts[sample][proportion]))
    
    return counts


def get_read_counts (read_counts_fn, samples):
    
    read_counts = defaultdict(lambda: defaultdict(int))
    read_counts_fn = os.path.realpath(read_counts_fn)
    
    err_suff = "\nSkipping plot with read counts..\n"
    
    # Check that the logs directory exists
    if not os.path.isfile(read_counts_fn):
        logging.error("Error - can't read read_counts file: {}{}".format(log_dir, err_suff))
        return None
    
    # Read the counts file
    try:
        with open(read_counts_fn, 'r') as fh:
            for line in fh:
                [fn, reads] = line.split("\t", 1)
                fn_base = os.path.basename(fn)
                fn_name = os.path.splitext(fn_base)[0]
                match = parse_filename(fn_name)
                if match is not None:
                    read_counts[match['sample']][match['proportion']] = reads
                else:
                    # We check we have everything at the end, don't panic!
                    continue
    except IOError as e:
        logging.error("Error - could not find cufflinks log file: {}\n{}{}".format(log_path, e, err_suff))
        return None
    
    # Loop through the sample names and check we have the reads
    for sample in samples:
        for proportion in samples[sample]:
            try:
                read_counts[sample][proportion]
            except NameError:
                logging.error("Error - couldn't find reads for {}, {}{}".format(sample, proportion, err_suff))
                return None
            
    return read_counts



def plot_gene_counts (counts, fpkm_cutoff, read_counts=None, output_fn='gene_counts'):
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
    min_x = 999999999
    max_x = 0
    for sample in sorted(counts):
        
        x_axis_values = []
        gene_counts = []
        
        for proportion in sorted(counts[sample]):
            gene_counts.append(counts[sample][proportion])
            if read_counts is not None:
                x_axis_values.append(int(read_counts[sample][proportion])/1000000)
            else:
                x_axis_values.append(float(proportion))
            min_x = min(min_x, min(x_axis_values))
            max_x = max(max_x, max(x_axis_values))
        axes.plot(x_axis_values, gene_counts, label=sample, color=colours[i], marker="x", markersize=3)
        i += 1
    
    # Tidy axes
    axes.set_xlim(min_x, max_x)
    axes.grid(True, zorder=0, which='both', axis='y', linestyle='-', color='#EDEDED', linewidth=1)
    axes.set_axisbelow(True)
    axes.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)
    
    # Make x axis proportions a percentage scale
    if read_counts is None:
        axes.set_xticklabels(["%d%%" % (z*100) for z in axes.get_xticks()])
    
    # Labels
    if read_counts is None:
        plt.xlabel('Proportion of sample')
    else:
        matplotlib.rcParams['mathtext.default'] = 'regular'
        plt.xlabel(r'Subsampled Read Counts ($\times 10^6$)')
    plt.ylabel("Number of genes with FPKM > {}".format(fpkm_cutoff))
    plt.title('Subsampled Gene Observations')

    # Legend
    axes.legend(loc='upper left', bbox_to_anchor = (1.02, 1.02), fontsize=8, markerscale=0)
    
    # SAVE OUTPUT
    png_fn = "{}.png".format(output_fn)
    pdf_fn = "{}.pdf".format(output_fn)
    logging.info("Saving to {} and {}".format(png_fn, pdf_fn))
    plt.savefig(png_fn)
    plt.savefig(pdf_fn)

    # Return the filenames
    return {'png': png_fn, 'pdf': pdf_fn}


    
    
if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("Plot number of observed genes at increasing read depths")
    parser.add_argument("-f", "--fpkm-cutoff", dest="fpkm_cutoff", default=0,
                        help="Cutoff at which to count genes as observed. Default: 0")
    parser.add_argument("-c", "--read-counts", dest="read_counts_fn", default=None,
                        help="File with read counts, used for x axis. See README for further info.")
    parser.add_argument("-o", "--output", dest="output_fn", default='gene_counts',
                        help="Plot output filename base. Default: gene_counts.png / .pdf")
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