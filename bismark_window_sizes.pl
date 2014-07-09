#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use POSIX qw(strftime);
use FindBin qw($Bin);
use lib "$FindBin::Bin/source";

my @bin_sizes = (100,200,300,400,500,1000,1500,2000,3000,4000,5000,10000,20000,30000,40000,50000,100000,200000,300000,400000,500000,1000000,2000000);
my @min_counts = (1,2,3,4,5,10);

my $only_strand_string;
my $only_strand;
my $regions;
my $min_coverage = 10;
my $min_counts_string;
my $bin_sizes_string;
my $append = "_windowSizes.txt";
my $quiet;
my $help;
my $config_result = GetOptions(
	"strand=s" => \$only_strand_string,
	"regions=s" => \$regions,
	"coverage=i" => \$min_coverage,
	"min_counts=s" => \$min_counts_string,
	"window_sizes=s" => \$bin_sizes_string,
	"append=s" => \$append,
	"quiet" => \$quiet,
	"help" => \$help
);
if($help){
	print "\n\n\tUSAGE: bismark_window_sizes.pl <coverage_file.cov>\n\n
--strand <for/rev> (default = both)
	Consider CpGs on only one strand - for (+) or rev (-)

--regions <regions.bed>
	BED file containing regions of interest

--coverage (default = 10)
	Coverage threshold to use for a window to be counted

--min_counts (default = 1,2,3,4,5,10)
	Comma separated list of minimum counts: the number of cytosines
	which pass the minimum coverage for a window to be counted

--window_sizes (default = 100,200,300,400,500,1000,1500,2000,3000,4000,5000,10000,20000,30000,40000,50000,100000,200000,300000,400000,500000,1000000,2000000)
	Comma separated list of window sizes to consider (in bp)
	
--append (default = _windowSizes.txt)
	Number of lines to process
	
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

if($only_strand_string and $only_strand_string ne 'for' and $only_strand_string ne 'rev'){
	die "Error! Strand must be either 'for' or 'rev' : $only_strand_string\n";
} elsif($only_strand_string){
	if($only_strand_string eq 'for'){
		$only_strand = '+';
	}
	if($only_strand_string eq 'rev'){
		$only_strand = '-';
	}
}

$min_coverage =~ s/\D+//g;
if(length($min_coverage) == 0 or $min_coverage < 1){
	die "Error! Coverage must be at least 1: $min_coverage\n";
}

if($min_counts_string){
	$min_counts_string =~ s/[\D,]//g;
	@min_counts = split(',', $min_counts_string);
	if(scalar @min_counts < 1){
		die "Error! Must have at least one minimum count: ".join(", ",@min_counts)."\n";
	}
}

if($bin_sizes_string){
	$bin_sizes_string =~ s/[\D,]//g;
	@bin_sizes = split(',', $bin_sizes_string);
	if(scalar @bin_sizes < 1){
		die "Error! Must have at least one window size: ".join(", ",@bin_sizes)."\n";
	}
}

my $min_bin_size = 99999999999999999;
my $max_bin_size = 0;
foreach my $binsize (@bin_sizes){
	if($min_bin_size > $binsize){
		$min_bin_size = $binsize;
	}
	if($max_bin_size < $binsize){
		$max_bin_size = $binsize;
	}
}

if(!$quiet) {
	warn "Configuration Options:\n";
	if($regions){ warn "\tROI File: $regions (Ignoring all reads outside of these regions)\n"; }
	if($only_strand){ warn "\tUsing Strand: $only_strand_string ($only_strand) (Ignoring all reads on other strands)\n"; }
	warn "\tCoverage Threshold: $min_coverage\n\tMinimum Counts: ".join(", ", @min_counts)."\n\tWindow Sizes: ".join(", ", @bin_sizes)."\n\tOutput filename append: $append\n\n";
}

my @filenames = @ARGV;
if(scalar @filenames == 0){ die "Error! No input filenames found.. Use --help for instructions\n"; }

# Print started
warn "Starting analysis at ".strftime("%H:%M:%S %a %b %e %Y", localtime)."\n";
my $started = time();

# Load regions of interest into memory
my %roi;
my %roi_windows;
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
		my $mb = int($start / 1000000);
		$roi{$chr}{$mb}{$start} = $end;
		
		# Counts for bins
		foreach my $binsize (@bin_sizes){
			my $bin_start = int($start / $binsize);
			my $bin_end = int($end / $binsize);
			# Make sure that start is before end
			if($bin_end < $bin_start){
				my $tmp = $bin_start;
				$bin_start = $bin_end;
				$bin_end = $tmp;
			}
			# Go through covered windows setting hash keys
			for(my $j = $bin_start; $j <= $bin_end; $j++){
				$roi_windows{$binsize}{$chr}{$j} = 1;
			}
		}
	}
	close ROI;
	warn "  ..loaded\n\n" unless($quiet);
}

# Go through each file
my $chr_col = 0;
my $pos_col = 1;
my $strand_col = 2;
my $methcount_col = 3;
my $unmethcount_col = 4;
foreach my $fn (@filenames){
	warn "Processing $fn\n" unless($quiet);
	open(IN, '<', $fn) or die "Couldn't open input file $fn: $? $!\n";
	
	# Go through each cytosine
	my $i = 0;
	my %bin_maxcoords;
	my %bincounts;
	while(<IN>){
		$i++;
		if($i % 1000000 == 0){
			warn "  ..$i lines\n" unless($quiet);
		}
		chomp;
		my @sections = split(/\t/);
		# What format input do we have?
		if($i == 1){
			# Genome wide cytosine report
			if($sections[2] eq '+' or $sections[2] eq '-'){
				# defaults already apply
			# bedGraph coverage file
			} elsif($sections[2] =~ /^\d+$/ && scalar @sections == 6){
				if($only_strand){
					die "Error! Stranded count specified but bedGraph coverage files don't specify strand...\n";
				}
				$chr_col = 0;
				$pos_col = 1;
				$strand_col = 0; # Not used
				$methcount_col = 4;
				$unmethcount_col = 5;
			}
		}
		my $chr = $sections[$chr_col];
		my $pos = $sections[$pos_col];
		my $strand = $sections[$strand_col];
		my $meth = $sections[$methcount_col];
		my $unmeth = $sections[$unmethcount_col];
		my $numcalls = $meth + $unmeth;
		my $mb = int($pos / 1000000);
		$chr =~ s/chr//;
		
		# Skip if the strand is wrong
		next if ($only_strand and $strand ne $only_strand);
		
		# Skip if the coverage is below threshold
		next if $numcalls < $min_coverage;
		
		# If we're using ROI, skip if not within region
		if($regions){
			my $covered = 0;
			foreach my $r_start ( keys(%{$roi{$chr}{$mb}})) {
				if($r_start < $pos and $roi{$chr}{$mb}{$r_start} > $pos){
					$covered = 1;
					last;
				} 
			}
			next unless $covered;
		}
		
		# Go through bins
		foreach my $binsize (@bin_sizes){
			my $bin_coord = int($pos / $binsize);
			$bincounts{$binsize}{$chr}{$bin_coord}++;
			if(!defined($bin_maxcoords{$binsize}{$chr}) or $pos > $bin_maxcoords{$binsize}{$chr}){
				$bin_maxcoords{$binsize}{$chr} = $pos;
			}
		}
		
	}
	close IN;
	
	# Work out the coverage bin stuff
	my @plot;
	my %percentages;
	my $binstats .= "Window Size\tTotal Windows";
	foreach my $count (@min_counts){
		$binstats .= "\tWindows Covered (min $count counts)\tPercentage Covered (min $count counts)";
	}
	$binstats .= "\n";
	foreach my $binsize (@bin_sizes){
		my $num_windows = 0;
		my %num_observed_windows;
		
		# COUNT TOTAL WINDOWS
		# If we have regions of interest
		if($regions){
			foreach my $chr (keys %{$roi_windows{$binsize}}) {
				$num_windows += scalar keys %{$roi_windows{$binsize}{$chr}};
			}
		# Else - use maximum observed co-ordinates
		} else {
			foreach my $chr (keys %{$bin_maxcoords{$binsize}}) {
				$num_windows += int($bin_maxcoords{$binsize}{$chr} / $binsize) + 1;
			}
		}
		
		# COUNT OBSERVED WINDOWS
		foreach my $chr (keys %{$bincounts{$binsize}}) {
			# Count windows with enough counts
			foreach my $cov_count (values %{$bincounts{$binsize}{$chr}}) {
				foreach my $count (@min_counts){
					if($cov_count >= $count){
						$num_observed_windows{$count}++;
					}
				}
				
			}
		}
		
		# Work out a percentage (safely)
		my %percentages;
		foreach my $count (@min_counts){
			if($num_windows == 0){
				$percentages{$count} = 0;
			} else {
				$percentages{$count} = sprintf("%.2f", (($num_observed_windows{$count} / $num_windows) * 100));
			}
		}
		
		# Save data for plot
		push @{$plot[0]}, log10($binsize);
		my $k = 1;
		foreach my $count (@min_counts){
			push @{$plot[$k]}, $percentages{$count};
			$k++;
		}
		
		# Add to string to print
		$binstats .= "$binsize\t$num_windows";
		foreach my $count (@min_counts){
			$binstats .= "\t".$num_observed_windows{$count}."\t".$percentages{$count}."%";
		}
		$binstats .= "\n";
	}
	
	# Print summary to STDERR if we're not being quiet
	warn $binstats unless($quiet);
	
	# Write everything to the log file
	my $output_fn = $fn.$append;
	open (OUT, '>', $output_fn) or die "Can't create coverage stats file $output_fn: $? $!\n";
	print OUT $binstats;
	close OUT;
	
	# Plot data
	# Check whether the module GD::Graph:lines and colour is installed
	eval{
		require GD::Graph::linespoints;
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
		my $graph = GD::Graph::linespoints->new(1000,600);
		add_colour(col1=>[27,158,119]);
		add_colour(col2=>[217,95,2]);
		add_colour(col3=>[117,112,179]);
		add_colour(col4=>[231,41,138]);
		add_colour(col5=>[102,166,30]);
		add_colour(col6=>[230,171,2]);
		add_colour(col7=>[166,118,29]);
		add_colour(col8=>[102,102,102]);
		$graph->set_title_font("$Bin/OpenSans-Regular.ttf", 16);
		$graph->set_legend_font("$Bin/OpenSans-Regular.ttf", 10);
		$graph->set_x_label_font("$Bin/OpenSans-Regular.ttf", 12);
		$graph->set_y_label_font("$Bin/OpenSans-Regular.ttf", 12);
		$graph->set_x_axis_font("$Bin/OpenSans-Regular.ttf", 10);
		$graph->set_y_axis_font("$Bin/OpenSans-Regular.ttf", 10);
	    $graph->set(
			 x_label              => 'Window Size',
			 y_label              => '% Windows',
			 title                => "Window Size Coverage - Min ".$min_coverage."x Coverage",
			 line_width           => 2,
			 markers			  => [4],
	 
			 x_min_value          => log10($min_bin_size),
			 x_max_value          => log10($max_bin_size),
			 x_tick_number		  => 10,
			 x_label_position     => 0.5,
			 x_number_format      => \&pow10,
	 
			 y_min_value          => 0,
			 y_max_value          => 100,
			 y_tick_number        => 10,
			 y_label_skip         => 2,
			 y_number_format      => '%d%%',
	 
			 bgclr                => 'white',
			 transparent          => 0,
			 legend_placement     => 'RT',
			 legend_spacing       => 6,
			 legend_marker_width  => 24,
			 legend_marker_height => 18,
			 dclrs                => [ qw(col1 col2 col3 col4 col5 col6 col7 col8) ],
		) or die $graph->error;
		if(scalar @min_counts > 1){
			my @legends;
			foreach my $count (@min_counts){
				push @legends, "Min $count Cs per window";
			}
			$graph->set_legend(@legends);
		}
		my $gd = $graph->plot(\@plot) or die $graph->error;
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
sub pow10 {
	my $n = shift;
	my $val = int(10**$n);
	my $printval = "$val bp";
	if($val >= 1000000){
		$printval = sprintf("%.1f mbp", ($val / 1000000));
	} elsif($val >= 1000){
		$printval = sprintf("%.1f kbp", ($val / 1000));
	}
	return $printval;
}