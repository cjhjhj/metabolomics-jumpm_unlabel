#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

## Initialization
my ($ms2Path) = @ARGV;
$ms2Path =~ s/\/$//;
my $outFile1 = $ms2Path . ".spectrum_matches";	## .spectrum_matches file containing all target and decoy structures with mscores for each feature
my $outFile2 = $ms2Path . ".structure"; ## .structure file containing a unique (target or decoy structure
open (OUT1, ">", $outFile1) or die "Cannot open $outFile1\n";
open (OUT2, ">", $outFile2) or die "Cannot open $outFile2\n";
my $header = "FeatureNo\tMass\tFormula\tName\tStructure\tInChi\tType\tAdduct\t" . 
			"NumMeasuredPeaks\tNumTheoreticalPeaks\tNumMatchedPeaks\tMscore\tIntensityRatio\n";
print OUT1 $header;
print OUT2 $header;
 
## Sort .score files according to the feature number
my @scoreFiles = glob("$ms2Path/*.score");
my @featureNums;
foreach my $scoreFile (@scoreFiles) {
	my ($featureNum) = $scoreFile =~ /f(\d+)\.MS2\.score/;
	push (@featureNums, $featureNum); 
}
my @index = sort {$featureNums[$a] <=> $featureNums[$b]} (0..$#featureNums);
@featureNums = @featureNums[@index];
@scoreFiles = @scoreFiles[@index];

## Read each .score file and write to output files
for (my $i = 0; $i < scalar(@scoreFiles); $i++) {
	open (SCORE, "<", $scoreFiles[$i]) or die "Cannot open $scoreFiles[$i]\n";	
	my $bestEntry;	## The best structure which has the highest "mscore"
	my $bestMscore = 0;	
	while (<SCORE>) {
		chomp ($_);
		my $line = $_;
		next if ($line =~ "Index");

		## Replace the first element of $line (i.e. contents of .score file)
		## with 'feature number' 
		my @elems = split(/\t/, $line);
		$elems[0] = $featureNums[$i];
		$line = join("\t", @elems);
		
		## Write to .spectrum_matches file
		print OUT1 $line, "\n";
		
		## Choose the entry with the highest mscore (i.e. best entry)
		my $mscore = (split(/\t/, $line))[-2];	## mscore is the second last element
		if ($mscore == -0) {
			$mscore = 0;
		}
		if (!defined ($bestEntry)) {
			$bestEntry = $line;
			$bestMscore = $mscore;
		} else {
			if ($mscore > $bestMscore) {
				$bestEntry = $line;
				$bestMscore = $mscore;
			}
		}
	}
	close (SCORE);
	## Sometimes, in spite of the existence of .out file, 
	## DB-matched entry in .out file may not have structure information (SMILES)
	## This entry in .out file needs to be skipped
	if (defined $bestEntry) {
		print OUT2 $bestEntry, "\n";
	}	
}
close (OUT1);
close (OUT2);
