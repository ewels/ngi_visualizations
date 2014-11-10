# Gene Body Coverage

A simple script to take the output from the 
[RSeQC](http://rseqc.sourceforge.net/#genebody-coverage-py)
`geneBody_coverage.py` script and plot this on single graph (png and pdf).

## Example output
![gene body coverage](../master/examples/geneBodyCoverage.png)

## Usage
This tool is a python script and should be run on the command line as follows:
```bash
python plot_gene_body_coverage.py *.geneBodyCoverage.txt
```

## Parameters

Command Line Flag | Description
----------------- | -------------------- | -----------
`<coverage_files>` | Required.<br>List of cufflinks results directories
`-o`, `--output` | Default: `gene_counts`<br>Plot output filename base. Default: `gene_counts.png` / `.pdf`
`-l`, `--log` | Default: `info`<br>Level of log messages to display. Can be `debug`, `info` or `warning`.
`-u`, `--log-output` | Default: `stdout`<br>Log output filename.
`-h`, `--help` | Display the help.

## Dependencies
The scripts processes output from
[RSeQC](http://rseqc.sourceforge.net/#genebody-coverage-py) and is written
in Python.

The following Python libraries are required:

* [matplotlib](http://matplotlib.org/)
* argparse
* logging
* os


