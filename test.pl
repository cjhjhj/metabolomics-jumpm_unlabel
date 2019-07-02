#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

#my @charges = (0, 0, 0, 1, 1, 1, 2, 2, 2);
#my %freq;
#foreach (@charges) {
#	$freq{$_}++;
#}
##my @tmp = sort {$freq{$b} <=> $freq{$a} || $a <=> $b} keys %freq;
#my $charge = (sort {$freq{$b} <=> $freq{$a} || $a <=> $b} keys %freq)[0];
#print Dumper($charge);

#my @array = ("abc\n", "bcd\n", "ccd\n");
#chomp(@array);
#print "@array\n";

#my ($n, $k, $r, $x) = (29, 12, 16, 8);
#my $p = fisherExact($n, $k, $r, $x);
#print "$p\n";
#
#sub fisherExact {
#	my ($n, $k, $r, $x) = @_;
#	my $p = 0;
#	for (my $i = 0; $i < $x; $i++) {
#		if (($k - $i) >= 0 && ($r - $i) >= 0 && (($n - $r) - ($k - $i)) >= 0) {			
#			$p += hypergeometric($n, $k, $r, $i);
#		} else {
#			next;
#		}
#	}
#	$p = 1 - $p;	## One-tailed p-value
#	return ($p);
#}
#
#sub hypergeometric {
#	my ($n, $k, $r, $x) = @_;
#	my $pr = exp(choose($r, $x) + choose($n - $r, $k - $x) - choose($n, $k));	
#	return ($pr);
#}
#
#sub choose {
#	my ($n, $k) = @_;
#	my ($result, $j) = (0, 1);
#	return 0 if $k > $n || $k < 0;
#	$k = ($n - $k) if ($n - $k) < $k;
#	while ($j <= $k) {
#		$result += log($n--);
#		$result -= log($j++);
#	}
#	return $result;
#}

#my @a = (1,3,3,6,7,8,9);
#my $m = median(@a);
#print "$m\n";
#
#
#sub median {
#        my (@data) = @_;
#        return unless @data;
#        return $data[0] unless scalar(@data) > 1;
#        @data= sort {$a <=> $b} @data;
#        return $data[$#data / 2] if @data&1;
#        my $mid= scalar(@data) /2;
#        return ($data[$mid - 1] + $data[$mid]) / 2;
#}
