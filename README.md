Visualizations
==============

A collection of next-gen sequencing visualisation scripts.

## Count Biotypes

This script takes an aligned BAM file and counts the overlaps with
different biotype flags within a GTF annotation file.

Annotation GTF (Gene Transfer Format) files can contain information about
coding sequences within the genome. In addition to specifying feature type
and position, GTF features can have associated annotation fields. One such
field is *biotype*, typically denoted by the `gene_type` or `biotype` flag
(see the [GENCODE format specifications of biotype flags](http://www.gencodegenes.org/gencode_biotypes.html).)

To get an overview of where reads from a next-generation sequencing library
are aligned within your reference genome, it can be interesting to annotate
overlaps with different types of features - for instance `rRNA` genes,
`protein_coding` genes and `miRNA` transcripts. This script does just that,
generating plots which show the frequency with which different biotype labels
are overlapped and how these overlaps are distributed throughout different
alignment lengths.

The script is written in Python and can be run on the command line or imported into another python script.

### Usage

On the command line:
```bash
python count_biotypes.py -g <annotation.gtf> <aligned_1.bam> .. <aligned_n.bam>
```

Within a python script:
```python
import count_biotypes
count_biotypes.main(annotation_file_path_, input_bam_file_paths):
```

If importing, individual functions can be called for a more 
fine-grained approach:
```python
(ftrs, bt_cts, bt_lnths) = count_biotypes.parse_gtf_biotypes(annotation_file_path)
(counts, lengths, output) = count_biotype_overlaps (aligned_bam, ftrs, bt_cts, bt_lnths)
(bargraph_png_fn, bargraph_pdf_fn) = plot_bars(counts, fn_basename)
(hist_png_fn, hist_pdf_fn) = plot_epic_histogram (lengths, fn_basename)
```

### Example output
The following plots were generated from a Total Small RNA run in Human cells,
accession [SRR1304304](http://www.ncbi.nlm.nih.gov/sra/?term=SRR1304304)

![Biotype overlaps](https://raw.githubusercontent.com/ewels/visualizations/master/examples/SRR1304304_trimmed_aligned_biotypeCounts.png)

![Biotype lengths](https://raw.githubusercontent.com/ewels/visualizations/master/examples/SRR1304304_trimmed_aligned_biotypeLengths.png)

![Biotype length percentages](https://raw.githubusercontent.com/ewels/visualizations/master/examples/SRR1304304_trimmed_aligned_biotypeLengthPercentages.png)

### Parameters

Command line flags and `count_biotypes.main()` arguments shown (in order received by `main()`).

def main(annotation_file, input_bam_list, biotype_flag='gene_type', feature_type='exon', num_lines=10000000, quiet=0):

Command Line Flag | main() argument name | Description
----------------- | -------------------- | -----------
`-g`, `--genome-feature-file` | `annotation_file` | Required. Path to annotation file.
`<input_bam_list>` | `input_bam_list` | Required. List of paths to aligned BAM files.
`-b`, `--biotype-flag` | `biotype_flag` | Default: `gene_type` (will also look for any flag containing `biotype`). Name of annotation flag to collect biotype label from.
`-t`, `--genome-feature` | `feature_type` | Default: `exon`. Type of feature to inspect within GTF file.
`-n`, `--num-lines` | `num_lines` | Default: 10 million. Number of lines to read from aligned BAM file.
`-q`, `--quiet` | `quiet` | Default: off. Prevents status messages being printed to stderr.
