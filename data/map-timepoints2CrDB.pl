#!/usr/bin/perl

# input: CrosstalkDB output csv containing all datasets from mouse tissue study (Tvardovskiy 2017), with header
# maps time points and number of the biological replicate to CrDB identifiers
# output: CrosstalkDB output csv with additional columns for time points and the biological replicate numbers, including header

use strict;
use warnings;
use Data::Dumper; $Data::Dumper::Sortkeys = 1;

# create list of tissues

my @tissues = ("brain", "heart", "liver", "kidney");

# create list of timepoints

my @timepoints = (("3 months") x 4, ("5 months") x 4, ("10 months") x 4, ("18 months") x 2, ("24 months") x 2);

# create list of replicates

my @reps = ((1..4) x 3, (1..2) x 2);

# create hash of CrDB ids and time points
# print hash content to file to pass to new data prep function

my %map;
my $tissue_counter = my $timerep_counter = 0;

my $cond_file = "sample_conditions.tsv";
open(my $cond_fh, '>', $cond_file) or die("Could not open $cond_file, exiting $!");
print $cond_fh "dataset id\ttissue\ttime point\tbiological replicate\n";

for (my $crdb_nr = 62; $crdb_nr <= 129; $crdb_nr++) {

	next if $crdb_nr == 64 || $crdb_nr == 65 || $crdb_nr == 66 || $crdb_nr == 67;

	my $tissue = $tissues[$tissue_counter];
	my $time = $timepoints[$timerep_counter];
	my $rep = $reps[$timerep_counter];

	my $leadzero = sprintf("%06d", $crdb_nr);
	my $crdb = "CrDB$leadzero";

	$map{$crdb} = "$tissue,$time,$rep";
	
	print $cond_fh "$crdb\t$tissue\t$time\t$rep\n";

	$timerep_counter++;

	if ($crdb_nr == 81 || $crdb_nr == 97 || $crdb_nr == 113) {

		$tissue_counter++;
		$timerep_counter = 0;

	}

}

close($cond_fh);

# read input data and add time point and biological replicate

my $incsv = $ARGV[0];
open(my $in_fh, '<', $incsv) or die("Could not open $incsv, exiting $!");

while (my $inline = <$in_fh>) {

	chomp($inline);
	my $out;	

	if ($inline =~ /accession number/) {

		$out = "$inline,\"timepoint\",\"biological replicate\"";

	} else {

		$inline =~ /"(CrDB\d+)"/;
		my $crdb = $1;

		my ($tissue, $time, $rep) = split(',', $map{$crdb});

		my $info = $map{$crdb};

		$out = "$inline,\"$time\",\"$rep\"";

	}

	print STDOUT "$out\n";

}

close($in_fh);
