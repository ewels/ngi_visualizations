#!/usr/bin/python
"""
count_biotypes.py
Takes aligned reads and uses HTSeq to count overlaps with different
biotype feature flags. Plots biotype proportions and a histogram
of read lengths broken up by biotype.
Input: GTF annotation file and aligned BAM files
"""

from __future__ import print_function

import argparse
from collections import defaultdict
import HTSeq
import logging
import numpy
import os

# Import matplot lib but avoid default X environment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def count_biotypes(annotation_file, input_bam_list, biotype_flag='gene_type', feature_type='exon', num_lines=10000000):
    """
    Count the biotypes
    """
    # Sanity check - make sure input files exist
    if annotation_file:
        if not os.path.isfile(annotation_file):
            raise IOError("Fatal error - can't find annotation file {}".format(annotation_file))
    else:
        raise ValueError("Fatal error - annotation file not specified")
    for fname in input_bam_list:
        if not os.path.isfile(fname):
            raise IOError("Fatal error - can't find input file {}".format(fname))
    
    # Parse the GTF file
    biotype_annotation = parse_gtf_biotypes(annotation_file, biotype_flag, feature_type)
    
    # Process files
    for fname in input_bam_list:
        logging.info("Processing {}".format(fname))
        
        # Make a copy of the biotype dicts
        biotype_count_dict = biotype_annotation['biotype_count_dict'].copy()
        
        # Generate counts
        biotype_count_dict = count_biotype_overlaps(fname, biotype_annotation['selected_features'], biotype_count_dict, num_lines)
        
        # Plot bar graph
        plot_basename = os.path.splitext(os.path.basename(fname))[0]
        plot_title = "{} Biotype Alignments".format(feature_type.title())
        bargraph_fns = plot_bars(biotype_count_dict['biotype_counts'], plot_basename, plot_title, True)
        log_bargraph_fns = plot_bars(biotype_count_dict['biotype_counts'], plot_basename, plot_title, False)
        
        # Plot epic histogram
        plot_title = "Read Lengths Overlapping {}s".format(feature_type.title())
        hist_fns = plot_epic_histogram (biotype_count_dict['biotype_lengths'], plot_basename, plot_title, False)
        percent_hist_fns = plot_epic_histogram (biotype_count_dict['biotype_lengths'], plot_basename, plot_title, True)



def parse_gtf_biotypes(annotation_file, biotype_label='gene_type', count_feature_type='exon'):
    """
    Custom function that uses HTSeq core to analyse overlaps
    with annotation features of different biotypes.
    Returns a dict of biotypes with their counts, a dict
    of biotypes with lists of read lengths, the total number
    of aligned reads and a summary string ready for printing
    """
    # Set up filenames & objects
    annotation_file = os.path.realpath(annotation_file)
    gtffile = HTSeq.GFF_Reader( annotation_file )
            
    # Go through annotation
    # Help from http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html#tour
    logging.info("\nParsing annotation file {}".format(annotation_file))
    selected_features = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    ignored_features = 0
    used_features = 0
    biotype_counts = {}
    biotype_lengths = {}
    biotype_counts['no_overlap'] = 0
    biotype_counts['multiple_features'] = 0
    biotype_lengths['no_overlap'] = defaultdict(int)
    biotype_lengths['multiple_features'] = defaultdict(int)
    feature_type_counts = defaultdict(int)
    feature_type_biotype_counts = defaultdict(lambda: defaultdict(int))
    for i, feature in enumerate(gtffile):
        if i % 100000 == 0 and i > 0:
            logging.debug("{} lines processed..".format(i))
        # Collect features and initialise biotype count objects
        if feature.type == count_feature_type:
            # See if we have another annotation that sounds like biotype
            # eg. Human ensembl calls it gene_biotype
            if biotype_label not in feature.attr and biotype_label == 'gene_type':
                for attr in feature.attr:
                    if 'biotype' in attr:
                        logging.warning("\nChanging biotype label from {} to {}".format(biotype_label, attr))
                        biotype_label = attr
            
            # Initiate count object and add feature to selected_features set  
            if biotype_label in feature.attr:
                used_features += 1
                selected_features[ feature.iv ] += feature.attr[biotype_label]
                biotype_counts[ feature.attr[biotype_label] ] = 0
                biotype_lengths[ feature.attr[biotype_label] ] = defaultdict(int)
            else:
                ignored_features += 1
                
        # Collect general annotation stats
        feature_type_counts[feature.type] += 1
        if biotype_label in feature.attr:
            feature_type_biotype_counts[feature.type][feature.attr[biotype_label]] += 1
    
    logging.info("\n\n{} features with biotype: {}".format(count_feature_type, used_features))
    logging.info("{} features without biotype: {}".format(count_feature_type, ignored_features))
    logging.info("{} biotypes to be counted: {}".format(count_feature_type, ', '.join(biotype_counts.keys())))
    
    logging.info("\nBiotype stats found for all feature types (using attribute '{}'):".format(biotype_label))
    for ft in sorted(feature_type_biotype_counts.keys()):
        num_ft_bts = len(feature_type_biotype_counts[ft].keys())
        num_features = 0
        for c,d in feature_type_biotype_counts[ft].iteritems():
            num_features += d
        logging.info("    {:20}\t{:4} biotypes\t{:6} labelled features".format(ft, num_ft_bts, num_features))
    
    if(used_features == 0):
        raise ValueError('No features have biotypes!')
    
    return {'selected_features': selected_features, 'biotype_count_dict': {'biotype_counts': biotype_counts, 'biotype_lengths':biotype_lengths}}


def count_biotype_overlaps(aligned_bam, selected_features, biotype_count_dict, number_lines=10000000):
    """
    Go thorough an aligned bam, counting overlaps with biotype features
    """
   
    # Set up filenames & objects
    aligned_bam = os.path.realpath(aligned_bam)
    bamfile = HTSeq.BAM_Reader( aligned_bam )
    
    # Go through alignments, counting transcript biotypes
    logging.info("\nReading BAM file (will stop at {}): ".format(number_lines))
    aligned_reads = 0
    for i, alnmt in enumerate(bamfile):
        if i > int(number_lines):
            i -= 1
            logging.info("Reached {} lines in the aligned file, exiting..".format(number_lines))
            break
        if i % 1000000 == 0 and i > 0:
            logging.debug("{} lines processed..".format(i))
        
        if alnmt.aligned:
            aligned_reads += 1
            iset = None
            for iv2, step_set in selected_features[ alnmt.iv ].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.intersection_update( step_set )
            
            key = 'multiple_features'
            if len(iset) == 1:
                key = list(iset)[0]
            elif len(iset) == 0:
                key = 'no_overlap'
                    
            biotype_count_dict['biotype_counts'][key] += 1
            biotype_count_dict['biotype_lengths'][key][alnmt.iv.length] += 1
    
    logging.info ("\n{} overlaps found from {} aligned reads ({} reads total)" \
                    .format(aligned_reads-biotype_count_dict['biotype_counts']['no_overlap'], aligned_reads, i))
    logging.info ("{} reads had multiple feature overlaps\n" \
                    .format(biotype_count_dict['biotype_counts']['multiple_features']))
    
    
    # Make a string table out of the counts
    counts_string = 'Type\tRead Count\n'
    for biotype in sorted(biotype_count_dict['biotype_counts'], key=biotype_count_dict['biotype_counts'].get, reverse=True):
        if biotype_count_dict['biotype_counts'][biotype] == 0:
            continue
        counts_string += "{}\t{}{}".format(biotype, biotype_count_dict['biotype_counts'][biotype], os.linesep)
    
    # Save to file
    file_basename = os.path.splitext(os.path.basename(aligned_bam))[0]
    counts_file = "{}_biotypeCounts.txt".format(file_basename)
    try:
        with open(counts_file, 'w') as fh:
            print(counts_string, file=fh);
    except IOError as e:
        raise IOError(e)
    
    # Return the counts
    return biotype_count_dict





def plot_bars(biotype_counts, output_basename, title="Annotation Biotype Alignments", logx=True):
    """
    Plots bar graph of alignment biotypes using matplotlib pyplot
    Input: dict of biotype labels and associated counts
    Input: file basename
    Input: Plot title
    Input: logx (t/f)
    Returns filenames of PNG and PDF graphs
    """
    
    # SET UP VARIABLES
    bar_width = 0.8
    total_reads = 0
    plt_labels = []
    plt_values = []
    for biotype in sorted(biotype_counts, key=biotype_counts.get):
        if biotype_counts[biotype] == 0:
            continue
        total_reads += biotype_counts[biotype]
        if biotype == 'no_overlap':
            continue
        plt_labels.append(biotype)
        plt_values.append(biotype_counts[biotype])
        
    ypos = numpy.arange(1, len(plt_labels)+1)
    
    minx = 0
    if logx:
        minx = 1
    
    # SET UP OBJECTS
    fig = plt.figure()
    axes = fig.add_subplot(111)
    
    # PLOT GRAPH
    barlist = axes.bar(minx, bar_width, plt_values, ypos, log=logx, align='center', orientation='horizontal', linewidth=0) 
    
    # Give more room for the labels on the left and top
    plt.subplots_adjust(left=0.25,top=0.8, bottom=0.15)
    
    # MAKE SPECIAL CASES GREY
    if 'multiple_features' in plt_labels:
        case_index = plt_labels.index('multiple_features')
        barlist[case_index].set_color('#999999')
    
    # Y AXIS
    axes.set_yticks(ypos)
    axes.set_yticklabels(plt_labels)
    axes.tick_params(left=False, right=False)
    axes.set_ylim((0,len(plt_labels)+1))
    
    # X AXIS
    axes.grid(True, zorder=0, which='both', axis='x', linestyle='-', color='#DEDEDE', linewidth=1)
    axes.set_axisbelow(True)
    if logx:
        axes.set_xscale('log')
    
    # SECONDARY X AXIS
    ax2 = axes.twiny()
    if logx:
        ax2.set_xscale('log')
    ax2.set_xlim(axes.get_xlim())
    ax1_ticks = axes.get_xticks()
    # I have no idea why I have to get rid of these two elements....
    ax1_ticks = ax1_ticks[1:-1]
    ax2.set_xticks(ax1_ticks)
    ax2.set_xlabel('Percentage of Aligned Reads')
    
    # SECONDARY AXIS LABELS
    def percent_total(counts):
        y = [(x/total_reads)*100 for x in counts]
        return ["%.2f%%" % z for z in y]
    ax2_labels = percent_total(ax2.get_xticks())
    ax2.set_xticklabels(ax2_labels)    
    
    # LABELS
    axes.set_xlabel('Number of Reads')
    axes.set_ylabel('Biotype')
    plt.text(0.5, 1.2, title, horizontalalignment='center',
                fontsize=16, weight='bold', transform=axes.transAxes)
    plt.text(0.5, 1.15, output_basename, horizontalalignment='center',
                fontsize=10, weight='light', transform = axes.transAxes)
    axes.tick_params(axis='both', labelsize=8)
    ax2.tick_params(axis='both', labelsize=8)
    if 'no_overlap' in biotype_counts:
        no_overlap_string = "{} reads had no feature overlap ({:.1%} of all {} aligned reads)" \
                            .format(biotype_counts['no_overlap'],
                            # ensure that these are being treated as floats not ints
                            (biotype_counts['no_overlap'] + 0.0) / (total_reads + 0.0) \
                            , total_reads)
        plt.text(0.5, -0.2, no_overlap_string, horizontalalignment='center',
                    fontsize=8, transform = axes.transAxes)
    
    # SAVE OUTPUT
    if logx:
        png_fn = "{}_biotypeCounts_log.png".format(output_basename)
        pdf_fn = "{}_biotypeCounts_log.pdf".format(output_basename)
    else:
        png_fn = "{}_biotypeCounts.png".format(output_basename)
        pdf_fn = "{}_biotypeCounts.pdf".format(output_basename)
    
    logging.info("Saving to {} and {}".format(png_fn, pdf_fn))
    fig.savefig(png_fn)
    fig.savefig(pdf_fn)
    
    # Return the filenames
    return {'png': png_fn, 'pdf': pdf_fn}


def plot_epic_histogram(biotype_lengths, output_basename, title="Annotation Biotype Lengths", percentage=False):
    """
    Plot awesome histogram of read lengths, with bars broken up by feature
    biotype overlap
    Input: dict of biotype labels with dict of lengths and counts
    Input: file basename
    Input: output fn
    Returns filenames of PNG and PDF graphs
    """
     
    # FIND MAX AND MIN LENGTHS, SET UP READ LENGTHS ARRAY
    min_length = 9999
    max_length = 0
    no_overlap_counts = 0
    total_reads = 0
    feature_reads = 0
    bp_counts = defaultdict(int)
    for bt in biotype_lengths:
        for x in biotype_lengths[bt]:
            total_reads += biotype_lengths[bt][x]
            if bt == 'no_overlap':
                no_overlap_counts += biotype_lengths[bt][x]
            else:
                feature_reads += biotype_lengths[bt][x]
                bp_counts[x] += biotype_lengths[bt][x]
                min_length = min(x, min_length)
                max_length = max(x, max_length)    
    
    # CUT OFF EXTREME READ LENGTHS
    # Trim off top and bottom 1% read lengths
    cum_count = 0
    first_percentile = (feature_reads + 0.0)/100.0
    nninth_percentile  = ((feature_reads + 0.0)/100.0)*99
    for x in xrange(min_length, max_length):
        cum_count += bp_counts[x]
        if (cum_count + 0.0) < first_percentile:
            min_length = x
        if (cum_count + 0.0) > nninth_percentile:
            if max_length+0.0 > x+0.0:
                max_length = x        
    
    # SET UP PLOT
    fig = plt.figure()
    axes = fig.add_subplot(111)
    plt.subplots_adjust(top=0.8, bottom=0.15, right=0.7)
    x_ind = range(min_length, max_length)
    bar_width = 0.8
    
    # PREPARE DATA
    bars = {}
    for bt in biotype_lengths:
        # Skip reads with no overlap
        if bt == 'no_overlap':
            continue
        values = []
        bt_count = 0
        for x in xrange(min_length, max_length):
            if x not in bp_counts:
                bp_counts[x] = 0
            if x in biotype_lengths[bt]:
                bt_count += biotype_lengths[bt][x]
                values.append(biotype_lengths[bt][x])
            else:
                values.append(0)
        if bt_count == 0:
            continue
        bars[bt_count] = (bt, values)
    
    # PLOT BARS
    pt = {}
    i = 0
    last_values = [0]*(max_length - min_length)
    legend_labels = []
    cols = distinguishable_colours(len(bars))
    for (count, bar) in sorted(bars.items(), reverse=True):
        (bt, values) = bar
        if(percentage):
            for (key,var) in enumerate(values):
                bp_count = bp_counts[key+min_length]
                if bp_count == 0:
                    values[key] = 0
                else:
                    values[key] = ((var+0.0)/(bp_counts[key+min_length]+0.0))*100
        pt[bt] = axes.bar(x_ind, values, width=bar_width, bottom=last_values, align='center', color=cols[i], linewidth=0)
        legend_labels.append(bt)
        last_values = [last_values+values for last_values,values in zip(last_values, values)]
        i += 1
    
    # TIDY UP AXES
    if percentage:
        axes.set_ylim((0,100))
    axes.set_xlim((min_length-1,max_length))
    axes.grid(True, zorder=0, which='both', axis='y', linestyle='-', color='#EDEDED', linewidth=1)
    axes.set_axisbelow(True)
    axes.tick_params(which='both', labelsize=8, direction='out', top=False, right=False)
    
    # LABELS
    axes.set_xlabel('Read Length (bp)')
    if(percentage):
        axes.set_ylabel('Percentage of Overlaps')
    else:
        axes.set_ylabel('Read Count')
    plt.text(0.5, 1.2, title, horizontalalignment='center',
                fontsize=16, weight='bold', transform=axes.transAxes)
    plt.text(0.5, 1.15, output_basename, horizontalalignment='center',
                fontsize=10, weight='light', transform = axes.transAxes)
    if 'no_overlap' in biotype_lengths:
        no_overlap_string = "{} reads had no feature overlap ({:.1%} of all {} aligned reads)" \
                            .format(no_overlap_counts, ((no_overlap_counts + 0.0) / (total_reads + 0.0)) \
                            , total_reads)
        plt.text(0.5, -0.2, no_overlap_string, horizontalalignment='center',
                    fontsize=8, transform = axes.transAxes)
    
    # LEGEND
    axes.legend(legend_labels, loc='upper left', bbox_to_anchor = (1.02, 1.02), fontsize=8)
    
    # SAVE OUTPUT
    png_fn = "{}_biotypeLengths.png".format(output_basename)
    pdf_fn = "{}_biotypeLengths.pdf".format(output_basename)
    if(percentage):
        png_fn = "{}_biotypeLengthPercentages.png".format(output_basename)
        pdf_fn = "{}_biotypeLengthPercentages.pdf".format(output_basename)
    
    logging.info("Saving to {} and {}".format(png_fn, pdf_fn))
    fig.savefig(png_fn)
    fig.savefig(pdf_fn)
    
    # Return the filenames
    return {'png': png_fn, 'pdf': pdf_fn}




# Stolen from @wefer
# https://github.com/wefer/color_picker/blob/master/dist_colors.py
def distinguishable_colours(n_colours):
    defaults = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928']*10
    # Import the packages that we need
    try:
        from colormath.color_objects import sRGBcolor, Labcolor
        from colormath.color_conversions import convert_color
        from matplotlib.colors import rgb2hex
    # fall back on a simple list
    except ImportError as e:
        logging.warning("Couldn't import modules for colours: {}\nFalling back to defaults.".format(e))
        return defaults
    
    bg = [1,1,1]    #Assumes background is white

    #Generate 30^3 RGB triples to choose from.
    n_grid = 30    
    x = numpy.linspace(0,1,n_grid)
    R = numpy.array([x]*900).flatten()
    G = numpy.array([[i]*30 for i in x]*30).flatten()
    B = numpy.array([[i]*900 for i in x]).flatten()

    rgb = numpy.array([R,G,B]).T #27000 by 3 matrix 

    if n_colours > len(rgb)/3:    #>27000
        logging.warning("You can't distinguish that many colours, dingus. " + \
                        "Tried to import {} colours".format(n_colours))
        return defaults

    #Convert to Lab colourspace
    lab = numpy.array([list(convert_color(sRGBcolor(i[0], i[1], i[2]), Labcolor).get_value_tuple()) for i in rgb])
    bglab = list(convert_color(sRGBcolor(bg[0], bg[1], bg[2]), Labcolor).get_value_tuple())

    #Compute the distance from background to candicate colours
    arr_length = len(rgb)
    mindist2 = numpy.empty(arr_length)
    mindist2.fill(float('Inf'))
    dX = lab - bglab
    dist2 = numpy.sum(numpy.square(dX), axis=1)
    mindist2 = numpy.minimum(dist2, mindist2)    

    #Pick the colours
    colours = numpy.zeros((n_colours, 3))
    lastlab = bglab

    for i in range(n_colours):
        dX = lab - lastlab    #subtract last from all colours in the list
        dist2 = numpy.sum(numpy.square(dX), axis=1)
        mindist2 = numpy.minimum(dist2, mindist2)
        index = numpy.argmax(mindist2)
        colours[i] = rgb[index]
        lastlab = lab[index]

    hex_colours =  [rgb2hex(item) for item in colours]

    return hex_colours
    


if __name__ == "__main__":
    # Command line arguments
    parser = argparse.ArgumentParser("Count read overlaps with different biotypes.")
    parser.add_argument("-g", "--genome-feature-file", dest="annotation_file", required=True,
                        help="GTF/GFF genome feature file to use for annotation (must match reference file)")
    parser.add_argument("-t", "--genome-feature", dest="feature_type", default='exon',
                        help="Type of annotation feature to count")
    parser.add_argument("-b", "--biotype-flag", dest="biotype_flag", default='gene_type',
                        help="GTF biotype flag (default = gene_type or *biotype*)")
    parser.add_argument("-n", "--num-lines", dest="num_lines", default=10000000,
                        help="Number of alignments to query")
    parser.add_argument("-l", "--log", dest="log_level", default='info',
                        help="Logging level: debug / info / warning")
    parser.add_argument("input_bam_list", metavar='<BAM file>', nargs="+",
                        help="List of input BAM filenames")
    kwargs = vars(parser.parse_args())
    
    # Initialise logger
    numeric_log_level = getattr(logging, kwargs['log_level'].upper(), None)
    if not isinstance(numeric_log_level, int):
        raise ValueError("Invalid log level: {}".format(kwargs['log_level']))
    logging.basicConfig(filename='count_biotypes.log', format='', level=numeric_log_level)
    # Remove logging paramter
    kwargs.pop('log_level', None)
    
    # Call count_biotypes()
    count_biotypes(**kwargs)
    