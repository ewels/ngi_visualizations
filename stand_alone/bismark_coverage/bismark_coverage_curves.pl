#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use POSIX qw(strftime);
use FindBin qw($Bin);
use lib "$FindBin::Bin/source";

my %roi_coverage;
my %notroi_coverage;
my %coverage;

my $num_Cs = 0;
my $num_roi_Cs = 0;
my $num_notroi_Cs = 0;

my %TS_roi_coverage;
my %TS_notroi_coverage;
my %TS_coverage;
my $TS_num_Cs = 0;
my $TS_num_roi_Cs = 0;
my $TS_num_notroi_Cs = 0;
my %BS_roi_coverage;
my %BS_notroi_coverage;
my %BS_coverage;
my $BS_num_Cs = 0;
my $BS_num_roi_Cs = 0;
my $BS_num_notroi_Cs = 0;

my $stranded;
my $regions;
my $min_cov = 0;
my $max_cov;
my $binsize = 1;
my $numlines = 1000000;
my $append = "_coverageStats.txt";
my $quiet;
my $help;
my $config_result = GetOptions(
	"stranded" => \$stranded,
	"regions=s" => \$regions,
	"min_cov=i" => \$min_cov,
	"max_cov=i" => \$max_cov,
	"binsize=i" => \$binsize,
	"numlines=i" => \$numlines,
	"append=s" => \$append,
	"quiet" => \$quiet,
	"help" => \$help
);
if($help){
	print "\n\n\tUSAGE: bismark_coverage_curves.pl *.cov\n\n
--regions <regions.bed>
	BED file containing regions of interest

--stranded
	Split the report up into forward and reverse strands

--min_cov (default = 0)
	Minimum coverage to consider

--max_cov (default = 15, 50 with capture regions)
	Maximum coverage to consider

--binsize (default = 1)
	Bin size to use

--numlines (default = 1000000)
	Number of lines to process

--append (default = _coverageStats.txt)
	String to append to the output filenames

--quiet
	Do not print status messages

--help
	Print this help message

\n";
	exit;
}
if(!$config_result){
	die "Error! could not parse command line options..\n";
}

# Set the default max coverage
unless($max_cov){
	if($regions){
		$max_cov = 50;
	} else {
		$max_cov = 15;
	}
}

# Print the config options
if(!$quiet) {
	warn "Configuration Options:\n";
	if($regions){ warn "\tROI File: $regions\n"; }
	warn "\tMin Coverage: $min_cov\n\tMax Coverage: $max_cov\n\tBin Size: $binsize\n\tNumber lines to process: $numlines\n\tOutput filename append: $append\n\n";
}

# Load in the input files
my @filenames = @ARGV;
if(scalar @filenames == 0){ die "Error! No input filenames found.. Use --help for instructions\n"; }

# Print started
warn "Starting analysis at ".strftime("%H:%M:%S %a %b %e %Y", localtime)."\n";
my $started = time();

# Load regions of interest into memory
my %roi;
my $roi_regions = 0;
my $roi_bp = 0;
if($regions){
	warn "\nLoading regions of interest..\n" unless($quiet);
	open(ROI, '<', $regions) or die "Couldn't open regions of interest file $regions: $? $!\n";
	while(<ROI>){
		chomp;
		my @sections = split(/\t/);
		my ($chr, $start, $end) = ($sections[0], $sections[1], $sections[2]);
		$chr =~ s/chr//;
		$start =~ s/\D//g;
		$end =~ s/\D//g;
		if(length($chr) > 0 and length($start) > 0 and length($end) > 0){
			my $mb = int($start / 1000000);
			$roi{$chr}{$mb}{$start} = $end;
			$roi_regions++;
			$roi_bp = $end - $start;
		}
	}
	close ROI;
	warn "  ..loaded\n\n" unless($quiet);
	warn "Found $roi_regions capture regions covering $roi_bp bp..\n" unless($quiet);
}

# Go through each file
my $orig_regions = $regions;
foreach my $fn (@filenames){
	$regions = $orig_regions;
	warn "Processing $fn\n" unless($quiet);
	open(IN, '<', $fn) or die "Couldn't open input file $fn: $? $!\n";
	# Go through each cytosine
	my $i = 0;
	while(<IN>){
		$i++;
		if($i % 100000 == 0){
			warn "  ..$i lines\n" unless($quiet);
		}
		if($i > $numlines){
			$i--;
			warn "Reached $i lines. Skipping rest of coverage file..\n" unless($quiet);
			last;
		}
		chomp;
		my @sections = split(/\t/);
		my $chr = $sections[0];
		my $pos = $sections[1];
		my $strand = $sections[2];
		my $meth = $sections[3];
		my $unmeth = $sections[4];
		my $numcalls = $meth + $unmeth;
		my $mb = int($pos / 1000000);
		$chr =~ s/chr//;
		if($regions){
			my $covered = 0;
			foreach my $r_start ( keys(%{$roi{$chr}{$mb}})) {
				if($r_start < $pos and $roi{$chr}{$mb}{$r_start} > $pos){
					$covered = 1;
					last;
				}
			}
			if($covered){
				$roi_coverage{$numcalls}++;
				$num_roi_Cs++;
				if($strand eq '+'){
					$TS_roi_coverage{$numcalls}++;
					$TS_num_roi_Cs++;
				} elsif($strand eq '-'){
					$BS_roi_coverage{$numcalls}++;
					$BS_num_roi_Cs++;
				}
			} else {
				$notroi_coverage{$numcalls}++;
				$num_notroi_Cs++;
				if($strand eq '+'){
					$TS_notroi_coverage{$numcalls}++;
					$TS_num_notroi_Cs++;
				} elsif($strand eq '-'){
					$BS_notroi_coverage{$numcalls}++;
					$BS_num_notroi_Cs++;
				}
			}
		}
		$coverage{$numcalls}++;
		$num_Cs++;
		if($strand eq '+'){
			$TS_coverage{$numcalls}++;
			$TS_num_Cs++;
		} elsif($strand eq '-'){
			$BS_coverage{$numcalls}++;
			$BS_num_Cs++;
		}
	}
	close IN;

	# Print some stats
	warn "Found $num_Cs total Cyotosines. $TS_num_Cs + strand, $BS_num_Cs bottom strand.\n" unless($quiet);;
	if($regions){
		warn "$num_roi_Cs cyotosines were within region of interest, $num_notroi_Cs outside.\n" unless($quiet);;
	}

	# Check we had some covered reads
	if($regions and $num_roi_Cs == 0){
		warn "\nError! No captured Cytosines found. Ignoring --regions\n";
		$regions = 0;
	}

	# Work out the coverage bin stuff
	my $binstats .= "Coverage Bin\tCount";
	if($stranded){
		$binstats .= "\tCount (+)\tCount (-)";
	}
	if($regions){
		if($stranded){
			$binstats .= "\tCount within ROI\tCount within Regions (+)\tCount within Regions (-)\tCount outside ROI\tCount outside Regions (+)\tCount outside Regions (-)";
		} else {
			$binstats .= "\tCount within ROI\tCount outside ROI";
		}
	}
	$binstats .= "\n";
	my %cov_counts;
	my %region_cov_counts;
	my %region_notcov_counts;
	my $binned_cov_counts = 0;
	my $binned_region_cov_counts = 0;
	my $binned_region_notcov_counts = 0;
	my %TS_cov_counts;
	my %TS_region_cov_counts;
	my %TS_region_notcov_counts;
	my $TS_binned_cov_counts = 0;
	my $TS_binned_region_cov_counts = 0;
	my $TS_binned_region_notcov_counts = 0;
	my %BS_cov_counts;
	my %BS_region_cov_counts;
	my %BS_region_notcov_counts;
	my $BS_binned_cov_counts = 0;
	my $BS_binned_region_cov_counts = 0;
	my $BS_binned_region_notcov_counts = 0;
	for (my $bin = $min_cov; $bin < $max_cov; $bin += $binsize){
		$cov_counts{$bin} = 0;
		$region_cov_counts{$bin} = 0;
		$region_notcov_counts{$bin} = 0;
		$TS_cov_counts{$bin} = 0;
		$TS_region_cov_counts{$bin} = 0;
		$TS_region_notcov_counts{$bin} = 0;
		$BS_cov_counts{$bin} = 0;
		$BS_region_cov_counts{$bin} = 0;
		$BS_region_notcov_counts{$bin} = 0;
		for (my $cov = $bin; $cov < ($bin + $binsize); $cov++){
			if(defined($coverage{$cov})){
				$cov_counts{$bin} += $coverage{$cov};
				$binned_cov_counts += $coverage{$cov};
			}
			if(defined($TS_coverage{$cov})){
				$TS_cov_counts{$bin} += $TS_coverage{$cov};
				$TS_binned_cov_counts += $TS_coverage{$cov};
			}
			if(defined($BS_coverage{$cov})){
				$BS_cov_counts{$bin} += $BS_coverage{$cov};
				$BS_binned_cov_counts += $BS_coverage{$cov};
			}
			if($regions){
				if(defined($roi_coverage{$cov})){
					$region_cov_counts{$bin} += $roi_coverage{$cov};
					$binned_region_cov_counts += $roi_coverage{$cov};
				}
				if(defined($TS_roi_coverage{$cov})){
					$TS_region_cov_counts{$bin} += $TS_roi_coverage{$cov};
					$TS_binned_region_cov_counts += $TS_roi_coverage{$cov};
				}
				if(defined($BS_roi_coverage{$cov})){
					$BS_region_cov_counts{$bin} += $BS_roi_coverage{$cov};
					$BS_binned_region_cov_counts += $BS_roi_coverage{$cov};
				}
				if(defined($notroi_coverage{$cov})){
					$region_notcov_counts{$bin} += $notroi_coverage{$cov};
					$binned_region_notcov_counts += $notroi_coverage{$cov};
				}
				if(defined($TS_notroi_coverage{$cov})){
					$TS_region_notcov_counts{$bin} += $TS_notroi_coverage{$cov};
					$TS_binned_region_notcov_counts += $TS_notroi_coverage{$cov};
				}
				if(defined($BS_notroi_coverage{$cov})){
					$BS_region_notcov_counts{$bin} += $BS_notroi_coverage{$cov};
					$BS_binned_region_notcov_counts += $BS_notroi_coverage{$cov};
				}
			}
		}
		my $line = join("\t", ($bin, $cov_counts{$bin}));
		if($stranded){
			$line .= join("\t", ($TS_cov_counts{$bin}, $BS_cov_counts{$bin}));
		}
		if($regions){
			if($stranded){
				$line .= "\t".join("\t", ($region_cov_counts{$bin}, $TS_region_cov_counts{$bin}, $BS_region_cov_counts{$bin},
										  $region_notcov_counts{$bin}, $TS_region_notcov_counts{$bin}, $BS_region_notcov_counts{$bin}));
			} else {
				$line .= "\t".join("\t", ($region_cov_counts{$bin}, $region_notcov_counts{$bin}));
			}
		}
		$binstats .= $line."\n";
	}

	# Sanity check for strandedness
	if($TS_num_Cs == 0 or $BS_num_Cs == 0){
		warn sprintf("\nError in total number of Top Strand Cs (%d) or Bottom Strand Cs (%d).\nCould be an input format problem? Ignoring --stranded..\n", $TS_num_Cs, $BS_num_Cs);
		$stranded = 0;
	}

	# Calculate summary & plot data
	my @plot;
	my @cov_counts_array;
	my @stranded_plot;
	my @TS_cov_counts_array;
	my @BS_cov_counts_array;
	my @regions_plot;
	my @region_cov_counts_array;
	my @region_notcov_counts_array;
	my @stranded_regions_plot;
	my @TS_region_cov_counts_array;
	my @TS_region_notcov_counts_array;
	my @BS_region_cov_counts_array;
	my @BS_region_notcov_counts_array;

	my $cum_count = $num_Cs;
	foreach my $bin (sort {$a <=> $b} keys(%cov_counts)){
		my $cov_percent = ($cum_count / $num_Cs)*100;
		push @{$plot[0]}, $bin;
		push @{$plot[1]}, $cov_percent;
		$cum_count -= $cov_counts{$bin};
		for(my $j = 0; $j < $cov_counts{$bin}; $j++){
			push @cov_counts_array, $bin;
		}
	}

	if($stranded){
		my $TS_cum_count = $TS_num_Cs;
		my $BS_cum_count = $BS_num_Cs;
		foreach my $bin (sort {$a <=> $b} keys(%TS_cov_counts)){
			my $TS_cov_percent = ($TS_cum_count / $num_Cs)*100;
			my $BS_cov_percent = ($BS_cum_count / $num_Cs)*100;
			push @{$stranded_plot[0]}, $bin;
			push @{$stranded_plot[1]}, $TS_cov_percent;
			push @{$stranded_plot[2]}, $BS_cov_percent;
			$TS_cum_count -= $TS_cov_counts{$bin};
			$BS_cum_count -= $BS_cov_counts{$bin};
			for(my $j = 0; $j < $TS_cov_counts{$bin}; $j++){
				push @TS_cov_counts_array, $bin;
			}
			for(my $j = 0; $j < $BS_cov_counts{$bin}; $j++){
				push @BS_cov_counts_array, $bin;
			}
		}
	}
	if($regions){
		my $cum_roi_count = $num_roi_Cs;
		my $cum_notroi_count = $num_notroi_Cs;
		foreach my $bin (sort {$a <=> $b} keys(%region_cov_counts)){
			my $cov_percent = ($cum_roi_count / $num_roi_Cs)*100;
			my $notcov_percent = ($cum_notroi_count / $num_notroi_Cs)*100;
			push @{$regions_plot[0]}, $bin;
			push @{$regions_plot[1]}, $cov_percent;
			push @{$regions_plot[2]}, $notcov_percent;
			$cum_roi_count -= $region_cov_counts{$bin};
			$cum_notroi_count -= $region_notcov_counts{$bin};
			for(my $j = 0; $j < $region_cov_counts{$bin}; $j++){
				push @region_cov_counts_array, $bin;
			}
			for(my $j = 0; $j < $region_notcov_counts{$bin}; $j++){
				push @region_notcov_counts_array, $bin;
			}
		}
	}
	if($stranded && $regions){
		my $TS_cum_roi_count = $TS_num_roi_Cs;
		my $TS_cum_notroi_count = $TS_num_notroi_Cs;
		my $BS_cum_roi_count = $BS_num_roi_Cs;
		my $BS_cum_notroi_count = $BS_num_notroi_Cs;
		foreach my $bin (sort {$a <=> $b} keys(%TS_region_cov_counts)){
			my $TS_cov_percent = ($TS_cum_roi_count / $TS_num_roi_Cs)*100;
			my $TS_notcov_percent = ($TS_cum_notroi_count / $TS_num_notroi_Cs)*100;
			my $BS_cov_percent = ($BS_cum_roi_count / $BS_num_roi_Cs)*100;
			my $BS_notcov_percent = ($BS_cum_notroi_count / $BS_num_notroi_Cs)*100;
			push @{$stranded_regions_plot[0]}, $bin;
			push @{$stranded_regions_plot[1]}, $TS_cov_percent;
			push @{$stranded_regions_plot[2]}, $BS_cov_percent;
			push @{$stranded_regions_plot[3]}, $TS_notcov_percent;
			push @{$stranded_regions_plot[4]}, $BS_notcov_percent;
			$TS_cum_roi_count -= $TS_region_cov_counts{$bin};
			$TS_cum_notroi_count -= $TS_region_notcov_counts{$bin};
			$BS_cum_roi_count -= $BS_region_cov_counts{$bin};
			$BS_cum_notroi_count -= $BS_region_notcov_counts{$bin};
			for(my $j = 0; $j < $TS_region_cov_counts{$bin}; $j++){
				push @TS_region_cov_counts_array, $bin;
			}
			for(my $j = 0; $j < $TS_region_notcov_counts{$bin}; $j++){
				push @TS_region_notcov_counts_array, $bin;
			}
			for(my $j = 0; $j < $BS_region_cov_counts{$bin}; $j++){
				push @BS_region_cov_counts_array, $bin;
			}
			for(my $j = 0; $j < $BS_region_notcov_counts{$bin}; $j++){
				push @BS_region_notcov_counts_array, $bin;
			}
		}
	}

	# Create summary
	my $summary = "\nSummary of coverage from $fn\n";
	$summary .= sprintf("\tCovered Cs: %d out of %d (%.2f%%)\n", $binned_cov_counts, $num_Cs, (($binned_cov_counts/$num_Cs)*100));
	if($regions){
		$summary .= sprintf("\tCovered Regions Cs: %d out of %d (%.2f%%)\n", $binned_region_cov_counts, $num_roi_Cs, (($binned_region_cov_counts/$num_roi_Cs)*100));
		$summary .= sprintf("\tCovered Cs outside Regions: %d out of %d (%.2f%%)\n", $binned_region_notcov_counts, $num_notroi_Cs, (($binned_region_notcov_counts/$num_notroi_Cs)*100));
	}
	$summary .= sprintf("\tMean coverage: %.2f\n\tMedian coverage: %d\n", mean(\@cov_counts_array), median(\@cov_counts_array));
	if($regions){
		$summary .= sprintf("\tMean Regions coverage: %.2f\n\tMedian Regions coverage: %d\n", mean(\@region_cov_counts_array), median(\@region_cov_counts_array));
		$summary .= sprintf("\tMean coverage outside Regions: %.2f\n\tMedian coverage outside Regions: %d\n", mean(\@region_notcov_counts_array), median(\@region_notcov_counts_array));
	}
	if($stranded){
		$summary .= "\nSummary of coverage on Top Strand (+):\n";
		$summary .= sprintf("\tCovered Cs: %d out of %d (%.2f%%)\n", $TS_binned_cov_counts, $TS_num_Cs, (($TS_binned_cov_counts/$TS_num_Cs)*100));
		if($regions){
			$summary .= sprintf("\tCovered Regions Cs: %d out of %d (%.2f%%)\n", $TS_binned_region_cov_counts, $TS_num_roi_Cs, (($TS_binned_region_cov_counts/$TS_num_roi_Cs)*100));
			$summary .= sprintf("\tCovered Cs outside Regions: %d out of %d (%.2f%%)\n", $TS_binned_region_notcov_counts, $TS_num_notroi_Cs, (($TS_binned_region_notcov_counts/$TS_num_notroi_Cs)*100));
		}
		$summary .= sprintf("\tMean coverage: %.2f\n\tMedian coverage: %d\n", mean(\@TS_cov_counts_array), median(\@TS_cov_counts_array));
		if($regions){
			$summary .= sprintf("\tMean Regions coverage: %.2f\n\tMedian Regions coverage: %d\n", mean(\@TS_region_cov_counts_array), median(\@TS_region_cov_counts_array));
			$summary .= sprintf("\tMean coverage outside Regions: %.2f\n\tMedian coverage outside Regions: %d\n", mean(\@TS_region_notcov_counts_array), median(\@TS_region_notcov_counts_array));
		}

		$summary .= "\nSummary of coverage on Bottom Strand (-):\n";
		$summary .= sprintf("\tCovered Cs: %d out of %d (%.2f%%)\n", $BS_binned_cov_counts, $BS_num_Cs, (($BS_binned_cov_counts/$BS_num_Cs)*100));
		if($regions){
			$summary .= sprintf("\tCovered Regions Cs: %d out of %d (%.2f%%)\n", $BS_binned_region_cov_counts, $BS_num_roi_Cs, (($BS_binned_region_cov_counts/$BS_num_roi_Cs)*100));
			$summary .= sprintf("\tCovered Cs outside Regions: %d out of %d (%.2f%%)\n", $BS_binned_region_notcov_counts, $BS_num_notroi_Cs, (($BS_binned_region_notcov_counts/$BS_num_notroi_Cs)*100));
		}
		$summary .= sprintf("\tMean coverage: %.2f\n\tMedian coverage: %d\n", mean(\@BS_cov_counts_array), median(\@BS_cov_counts_array));
		if($regions){
			$summary .= sprintf("\tMean Regions coverage: %.2f\n\tMedian Regions coverage: %d\n", mean(\@BS_region_cov_counts_array), median(\@BS_region_cov_counts_array));
			$summary .= sprintf("\tMean coverage outside Regions: %.2f\n\tMedian coverage outside Regions: %d\n", mean(\@BS_region_notcov_counts_array), median(\@BS_region_notcov_counts_array));
		}
	}


	$summary .= "\n";

	# Print summary to STDERR if we're not being quiet
	warn $summary unless($quiet);

	# Write everything to the log file
	my $output_fn = $fn.$append;
	open (OUT, '>', $output_fn) or die "Can't create coverage stats file $output_fn: $? $!\n";
	print OUT $summary."\n";
	print OUT $binstats;
	close OUT;


	# Check whether the module GD::Graph:lines and colour is installed
	eval{
		require GD::Graph::lines;
		require GD::Graph::colour;
		require GD::Image;
		GD::Graph::lines->import();
		GD::Graph::colour->import(qw(:colours :lists :files :convert));
		GD::Graph::Image->import();
	};
	if($@) { # something went wrong
		warn "Perl module GD::Graph::lines or GD::Graph::colour not found, skipping plots\n";
	} else {
		warn "Plotting coverage graph\n" unless ($quiet);
		my $graph_fn = $fn.$append.'.png';
		my $graph = GD::Graph::lines->new(800,600);
	    add_colour(blue=> [31,120,180]);
		add_colour(orange=> [255,127,0]);
	    add_colour(green=> [51,160,44]);
		add_colour(purple=> [106,61,154]);
		$graph->set_title_font("$Bin/OpenSans-Regular.ttf", 16);
		$graph->set_legend_font("$Bin/OpenSans-Regular.ttf", 12);
		$graph->set_x_label_font("$Bin/OpenSans-Regular.ttf", 12);
		$graph->set_y_label_font("$Bin/OpenSans-Regular.ttf", 12);
		$graph->set_x_axis_font("$Bin/OpenSans-Regular.ttf", 10);
		$graph->set_y_axis_font("$Bin/OpenSans-Regular.ttf", 10);
	    $graph->set(
			 x_label              => 'Minimum Coverage',
			 y_label              => '% CpGs',
			 title                => 'Coverage Decay Plot',
			 line_width           => 2,

			 x_min_value          => 0,
			 x_max_value          => $max_cov,
			 x_label_skip         => 5,
			 x_label_position     => 0.5,
			 x_number_format      => '%d x',

			 y_min_value          => 0,
			 y_max_value          => 100,
			 y_tick_number        => 10,
			 y_label_skip         => 2,
			 y_number_format      => '%d%%',

			 bgclr                => 'white',
			 transparent          => 0,
			 legend_placement     => 'BC',
			 legend_spacing       => 6,
			 legend_marker_width  => 24,
			 legend_marker_height => 18,
			 dclrs              => [ qw(purple green orange blue)],
		) or die $graph->error;
		my $gd;
		if(!$stranded && !$regions){
   		    $gd = $graph->plot(\@plot) or die $graph->error;
		} elsif($stranded && !$regions){
			$graph->set_legend('Top Strand (+)','Bottom Strand (-)');
		    $gd = $graph->plot(\@stranded_plot) or die $graph->error;
		} elsif(!$stranded && $regions){
			$graph->set_legend('Within ROI','Outside ROI');
		    $gd = $graph->plot(\@regions_plot) or die $graph->error;
		} elsif($stranded && $regions){
			$graph->set_legend('Within Regions (+)','Within Regions (-)','Outside Regions (+)','Outside Regions (-)');
		    $gd = $graph->plot(\@stranded_regions_plot) or die $graph->error;
		}
		# Save plot to file
	    open (PLOT,'>',$graph_fn) or die "Failed to write to file for plot $graph_fn: $!\n\n";
	    binmode PLOT;
	    print PLOT $gd->png;
		close PLOT;
	} # End of plotting

} # End of files loop

warn "Finished analysis at ".strftime("%H:%M:%S %a %b %e %Y", localtime)."\n" unless($quiet);
my $duration = sprintf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime (time() - $started))[7,2,1,0];
warn "Analysis took $duration\n" unless($quiet);


sub log10 {
	my $n = shift;
	return log($n)/log(10);
}
sub mean {
    my ($arrayref) = @_;
    my $result = 0;
    foreach (@$arrayref) {
		$result += $_;
	}
    return $result / @$arrayref;
}
sub median {
	$_[0]->[ @{$_[0]} / 2 ];
}
