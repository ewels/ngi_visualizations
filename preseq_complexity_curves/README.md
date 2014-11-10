
# preseq Complexity Curves
[Preseq](http://smithlabresearch.org/software/preseq/) is a piece of software
designed to estimate how library complexity varies with increasing sequencing
depth. It helps to decide whether further sequencing would yield further
information.

The output from preseq is purely text. This script takes a number of preseq
output files and plots them using the Python matplotlib package. Optionally,
a file can be provided with read counts, which will be plotted on the lines
to show to what depth the current libraries have been sequenced. This file
can contain two columns of numbers, for instance if unique reads have been
calculated using a different method (_eg._ Picard).

Note that preseq uses a different method to calculate how many unique molecules
there are to most duplicate removal tools. See [this seqanswers thread](http://seqanswers.com/forums/showthread.php?t=27798&page=2)
for more information.

## Example output
Plotting points with read counts versus interpolated unique read counts (from preseq):
![Just read counts](../examples/complexity_curves_readcounts.png)

Plotting points with read counts versus deduplicated read counts (using Picard):
![Read counts and unique reads](../examples/complexity_curves_PicardDups.png)

## Usage

This tool is a python script and should be run on the command line as follows:

```bash
python plot_complexity_curves.py *.preseq
```

### Plotting the "real" data
To the real data as points, a file contaning the read counts must be referenced
using the `-r` argument. You can generate such a file using [samtools](http://www.htslib.org/)
with the following command:
```
echo "Sample_1 "$(samtools view -c -F 4 Sample_1.bam))
```

If you would like to count two files - one with and one without duplicates,
you can execute the following:
```
echo "Sample_1 "$(samtools view -c -F 4 Sample_1.bam)" "$(samtools view -c -F 4 Sample_1_dedup.bam)
```

The file should have at least two columns - the first being the input filename
and the second being the read counts (space delimited). The third optional
column with unique counts is optional.

## Parameters

Arguments shown in order received by `plot_complexity_curves()`.

Command Line Flag | `plot_complexity_curves()` argument name | Description
----------------- | -------------------- | -----------
`<preseq_files>` | `ccurves` | Required.<br>List of paths to preseq output files.
`-r`, `--real-counts` | `real_counts_path` | Default: none.<br>Path to file with read counts.
`-o`, `--output-name` | `output_name` | Default: `complexity_curves`.<br>Output file name
`-m`, `--x-min` | `x_min` | Default: `0`.<br>Minimum x axis limit.
`-x`, `--x-max` | `x_max` | Default: 500 million.<br>Maximum x axis limit.
`-h`, `--help` | - | Display full help text.

## Dependencies

The script is written in Python. The following libraries are required:

* [matplotlib](http://matplotlib.org/)
* [numpy](http://www.numpy.org/)
* [pandas](http://pandas.pydata.org/)
* argparse
* glob
* os
* subprocess
* sys
* yaml


