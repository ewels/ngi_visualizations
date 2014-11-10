Visualizations
==============

A collection of next-gen sequencing visualisation scripts. Click a script's
name to go to it's subdirectory which will contain a detailed `README.md`
file with examples and instructions.

* [Count Biotypes](count_biotypes/)
	* Uses HTSeq to plot read overlaps with different feature biotype flags
* [preseq Complexity Curves](preseq_complexity_curves/)
* [Subsampled Gene Observations](subsampled_gene_observations/)
    * Group of scripts to plot the number of observed genes at varying sample
    subsampling proportions. Can give an impression of library complexity on
    a biological level.
* [Gene Body Coverage](gene_body_coverage/)
* [Bismark Addons](bismark/)
	* [Bismark Coverage Curves](bismark/#bismark-coverage-curves) - Plots the proportion of cytosines meeting increasing coverage thresholds
	* [Bismark Window Sizes](bismark/#bismark-window-sizes) - Plots the proportion of windows passing observation thresholds with increasing window sizes
* [Alignment Summaries](alignment_summaries/)
	* Two scripts to parse log files containing alignment stats from bowtie,
		bowtie 2 or tophat and generate overview HTML reports

See below for example outputs. Click an image to go to that script.


[<img width="45%" src="examples/SRR1304304_trimmed_aligned_biotypeCounts.png">](count_biotypes/)

[<img width="45%" src="examples/SRR1304304_trimmed_aligned_biotypeLengths.png">](count_biotypes/)

[<img width="45%" src="examples/complexity_curves_readcounts.png">](preseq_complexity_curves/)

[<img width="45%" src="examples/subsampled_gene_observations.png">](subsampled_gene_observations/)

[<img width="45%" src="examples/geneBodyCoverage.png">](gene_body_coverage/)

[<img width="45%" src="examples/coverageStats.png">](bismark/)

[<img width="45%" src="examples/windowSizes_roi.png">](bismark/)

[<img width="45%" src="examples/bowtie_align_screenshot.png">](alignment_summaries/)

