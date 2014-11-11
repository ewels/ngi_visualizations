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

### See below for example outputs. Click an image to go to that script.

<table>
  <tr>
    <th colspan="2"><a href="count_biotypes/">Count Biotypes</a></th>
  </tr>
  <tr>
    <td>
      <a href="count_biotypes/" title="Count Biotypes">
        <img src="examples/SRR1304304_trimmed_aligned_biotypeCounts.png">
      </a>
    </td>
    <td>
      <a href="count_biotypes/" title="Count Biotypes">
        <img src="examples/SRR1304304_trimmed_aligned_biotypeLengths.png">
      </a>
    </td>
  </tr>
</table>

<table>
  <tr>
    <th><a href="preseq_complexity_curves/">preseq Complexity Curves</a></th>
    <th><a href="subsampled_gene_observations/">Subsampled Gene Observations</a></th>
  </tr>
  <tr>
    <td>
      <a href="preseq_complexity_curves/" title="preseq Complexity Curves">
        <img src="examples/complexity_curves_readcounts.png">
      </a>
    </td>
    <td>
      <a href="subsampled_gene_observations/" title="Subsampled Gene Observations">
        <img src="examples/subsampled_gene_observations.png">
      </a>
    </td>
  </tr>
</table>

<table>
  <tr>
    <th><a href="gene_body_coverage/">Gene Body Coverage</a></th>
    <th><a href="alignment_summaries/">Alignment Summaries</a></th>
  </tr>
  <tr>
    <td>
      <a href="gene_body_coverage/" title="Gene Body Coverage">
        <img src="examples/geneBodyCoverage.png">
      </a>
    </td>
    <td>
      <a href="alignment_summaries/" title="Alignment Summaries">
        <img src="examples/bowtie_align_screenshot.png">
      </a>
    </td>
  </tr>
</table>

<table>
  <tr>
    <th><a href="bismark/#bismark-coverage-curves">Bismark Coverage Curves</a></th>
    <th><a href="bismark/#bismark-window-sizes">Bismark Window Sizes</a></th>
  </tr>
  <tr>
    <td>
      <a href="bismark/#bismark-coverage-curves" title="Bismark Coverage Curves">
        <img src="examples/coverageStats.png">
      </a>
    </td>
    <td>
      <a href="bismark/#bismark-window-sizes" title="Bismark Window Sizes">
        <img src="examples/windowSizes_roi.png">
      </a>
    </td>
  </tr>
</table>


