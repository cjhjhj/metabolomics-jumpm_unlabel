#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;

my $paramFile;
GetOptions('-p=s' => \$paramFile,);
if (!-e ($paramFile)) {
	print "Please input the parameter file\n\n";
	exit;
}
$paramFile = abs_path($paramFile);

## Load R module first
system ("module load R/3.5.1");
#my $cmd = "bsub -P prot -q standard -R \"rusage[mem=100000]\" -Ip _jumpm.pl -p " . join(" ", @ARGV);
my $cmd = "bsub -P prot -q standard -R \"rusage[mem=100000]\" -Ip testJumpm.pl -p " . join(" ", @ARGV);
system($cmd);
