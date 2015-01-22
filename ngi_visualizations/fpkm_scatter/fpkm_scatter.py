#!/usr/bin/python
"""
fpkm_scatter.py

Takes FPKM counts from two conditions and makes a scatter plot.
Also calculates a r-squared correlation score.

Either takes two Cufflinks FPKM summary files or a summary file
with multiple samples.
"""

from __future__ import print_function

import argparse
import logging
import numpy as np
import os
import re

# Import matplot lib but avoid default X environment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def make_fpkm_scatter_plots (input_files, summary=False, output_fn='gene_counts', linear=False):
    """
    Main function. Takes input files and makes a plot.
    """
    
    # Get proper paths for input files
    input_1 = os.path.realpath(input_files[0])
    input_2 = os.path.realpath(input_files[1])
    
    # First off, assume that we have two FPKM files from cufflinks
    if summary is False:
       
        # Parse the input files
        sample_1 = load_fpkm_counts(input_1)
        sample_2 = load_fpkm_counts(input_2)
        
        # What are the sample names?
        sample_1_name = os.path.splitext(os.path.basename(input_1))[0]
        sample_2_name = os.path.splitext(os.path.basename(input_2))[0]
        
        # File name
        if output_fn is None:
            output_fn = "{}-{}".format(sample_1_name, sample_2_name)

        # Make the plot
        plot_filenames = plot_fpkm_scatter(sample_1, sample_2, sample_1_name, sample_2_name, output_fn=output_fn, linear=linear)
    
    # We have a summary file
    else:
        # Parse the input files
        # returns counts[sample_name][gene_id] = fpkm_count
        condition_1 = load_summary_fpkm_counts(input_1)
        condition_2 = load_summary_fpkm_counts(input_2)
        
        # Condition basenames
        cond_1_basename = os.path.splitext(os.path.basename(input_1))[0]
        cond_2_basename = os.path.splitext(os.path.basename(input_2))[0]
        
        # Find project names from sample ID (P1234_*)
        pname_re = re.compile('^P\d+_')
        cond_1_proj = pname_re.search(condition_1.keys()[0]).group(0)
        cond_2_proj = pname_re.search(condition_2.keys()[0]).group(0)
        
        # Go through each sample pair
        for cond_1_sample in condition_1.keys():
            try:
                cond_2_sample = cond_1_sample.replace(cond_1_proj, cond_2_proj)
                condition_2[cond_2_sample]
                outfile = "{}_{}-{}_{}".format(cond_1_basename, cond_1_sample, cond_2_basename, cond_2_sample)
                
                plot_filenames = plot_fpkm_scatter(condition_1[cond_1_sample], condition_2[cond_2_sample], cond_1_sample, cond_2_sample, output_fn=outfile, linear=linear)
            except KeyError:
                logging.warning("Warning: Sample {} not found in {}".format(sample, cond_2_basename))
    
    
def load_fpkm_counts (file):
    """
    Reads FPKM file and returns a dict of counts.
    Input: Path to FPKM file
    Output: Dictionary. Keys are transcript names, values are counts.
    """
    
    counts = {}
    logging.info("Processing {}".format(file))
    if os.path.isfile(file):
        try:
            with open(file, 'r') as fh:
                # Skip the header
                next(fh)
                # Iterate through the file
                for line in fh:
                    line = line.strip()
                    cols = line.split("\t")
                    gene_id = cols[0]
                    FPKM = cols[9]
                    counts[gene_id] = FPKM
        except IOError as e:
            logging.error("Error loading FPKM file: {}".format(file))
            raise IOError(e)
    else:
        raise IOError("Couldn't find FPKM file: {}".format(file))
    
    return counts


def load_summary_fpkm_counts (file):
    """
    Reads FPKM summary file and returns counts.
    Input: Path to summary FPKM file
    Output: Nested Defaultdict. First keys are sample names,
            second keys are gene names, values are counts.
            counts[sample_name][gene_id] = fpkm_count
    """
    
    counts = {}
    
    # ENSEMBL_ID    Gene_ID    Sample_1    Sample_2    Sample_3    Sample_4
    # ENSG00000000003    TSPAN6    0.1234    0.1234    0.1234    0.1234
    
    # Read the counts file
    if os.path.isfile(file):
        try:
            with open(file, 'r') as fh:
                header = fh.readline()
                sample_names = header.split("\t")[2:]
                for name in sample_names:
                    counts[name] = {}
                for line in fh:
                    sections = line.split("\t")
                    gene_id = sections[0]
                    for i in range(2, len(sections)-2):
                        counts[sample_names[i]][gene_id] = sections[i]
        except IOError as e:
            logging.error("Error - could not find cufflinks log file: {}\n{}{}".format(log_path, e, err_suff))
            return None
    else:
        raise IOError("Couldn't find FPKM file: {}".format(file))
            
    return counts



def plot_fpkm_scatter (sample_1, sample_2, x_lab, y_lab, output_fn=False, linear=False):
    """
    Plots scatter plot of FPKM counts. Also calculates r-squared correlation
    and prints this to STDOUT as well as including it in the graph.
    Input: 2 x counts[gene_id] = fpkm_count and 2 x sample names
    Input: title=Plot title
    Returns a dict containing filenames of PNG and PDF graphs as a dict for 
    each sample.
    """
    
    # Output filename
    if output_fn is False:
        output_fn = "{}_{}".format(x_lab, y_lab)
    
    # SET UP PLOT
    fig = plt.figure()
    axes = fig.add_subplot(111, aspect=1.0)
    
    i = 0
    x_vals = []
    y_vals = []
    
    # Collect the paired FPKM counts
    missing_genes = 0
    for gene in sample_1.keys():
        try:
            sample_2[gene]
        except NameError:
            missing_genes += 1
        else:
            x_vals.append(float(sample_1[gene]))
            y_vals.append(float(sample_2[gene]))
    if missing_genes > 0:
        logging.warn("Warning: {} genes mentioned in sample 1 not found in sample 2".format(missing_genes))
    
    # Calculate the r squared    
    corr = np.corrcoef(x_vals, y_vals)[0,1]
    r_squared = corr ** 2
    print ("R squared for {} = {}".format(output_fn, r_squared))

    # Make the plot
    axes.plot(x_vals, y_vals, 'o', markersize=1)
    plt.xlabel(x_lab)
    plt.ylabel(y_lab)
    plt.title("FPKM Counts for {} and {}".format(x_lab, y_lab))
    
    # Axes scales
    x1,x2,y1,y2 = axes.axis()
    max_xy = max(x2, y2)
    if linear is True:
        axes.set(xscale='log', xlim=[0,max_xy])
        axes.set(yscale='log', ylim=[0,max_xy])        
    else:
        axes.set(xscale='log', xlim=[1,max_xy])
        axes.set(yscale='log', ylim=[1,max_xy])
    
    # Tidy axes
    axes.set_axisbelow(True)
    axes.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)
    
    # Add a label about r squared
    plt.subplots_adjust(bottom=0.15)
    plt.text(1, -0.15, r'$r^2$ = {:2f}'.format(r_squared),
                horizontalalignment='right', fontsize=8, transform=axes.transAxes)
    
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
    parser = argparse.ArgumentParser("Make a scatter plot of FPKM counts between conditions")
    parser.add_argument("-s", "--summary", dest="summary", action='store_true',
                        help="Input files are summary FPKM files.")
    parser.add_argument("-o", "--output", dest="output_fn", default=None,
                        help="Plot output filename base. Default: sample1_sample2.png / .pdf")
    parser.add_argument("-n", "--linear", dest="linear", action='store_true',
                        help="Plot using linear axes instead of log10")
    parser.add_argument("-l", "--log", dest="log_level", default='info', choices=['debug', 'info', 'warning'],
                        help="Level of log messages to display")
    parser.add_argument("-u", "--log-output", dest="log_output", default='stdout',
                        help="Log output filename. Default: stdout")
    parser.add_argument("input_files", metavar='<input files>', nargs=2,
                        help="List of summary FPKM / cufflinks FPKM results files. See README.MD for more information.")
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
    make_fpkm_scatter_plots(**kwargs)