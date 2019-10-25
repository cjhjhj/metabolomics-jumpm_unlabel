#!/usr/bin/perl

use strict;
use warnings;
use Cwd 'abs_path';
use Getopt::Long;
use FindBin qw($Bin);
use lib "$Bin";

my ($paramFile, $queue, $mem);
GetOptions('-p=s' => \$paramFile, '--queue=s'=>\$queue, '--mem=s'=>\$mem);
if (!-e ($paramFile)) {
	print "Please input the parameter file\n\n";
	exit;
}
$paramFile = abs_path($paramFile);
if(!defined($queue) && !defined($mem)) {
    $queue = 'standard';
    $mem = 10000;
}
elsif(!defined($queue) && defined($mem)) {
    print "\t--mem cannot be used without --queue\n";
    exit(1);
}
elsif(!defined($mem)) {
    $mem = 10000;
}
my $cmd = "bsub -P prot -q $queue -R \"rusage[mem=$mem]\" -Ip $Bin/_jumpm.pl -p $paramFile " . join(" ", @ARGV);
system($cmd);
