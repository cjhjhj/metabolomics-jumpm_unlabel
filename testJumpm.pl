#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

#my $nFiles = 4200;
#my $maxJobs = 200;
#my $filesPerJob;
#if ($nFiles <= $maxJobs) {
#	$maxJobs = $nFiles;
#	$filesPerJob = 1;
#} else {
#	$filesPerJob = int($nFiles / $maxJobs) + 1;
#}
#my $nTotalJobs = int($nFiles / $filesPerJob - 0.0001) + 1;
#my $nJobs = 0;
#for (my $i = 0; $i < $nTotalJobs; $i++) {
#	$nJobs++;		
#	for (my $j = 0; $j < $filesPerJob; $j++) {
#		my $k = $filesPerJob * $i + $j;
#		if ($k >= $nFiles) {
#			last;
#		}
#		print "job# $nJobs, file# $k\n";
#	}
#	print "$i\n";
#}

my $paramFile;
GetOptions('-p=s' => \$paramFile,);
if (!-e ($paramFile)) {
	print "Please input the parameter file\n\n";
	exit;
}
print "$paramFile\n";
print Dumper(@ARGV), "\n";
