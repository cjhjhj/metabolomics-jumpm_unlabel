#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $mzXML = "./IROA_IS_NEG_1/IROA_IS_NEG_1.mzXML";
my $tmpFeatureFile = "./IROA_IS_NEG_1/IROA_IS_NEG_1.tmp.feature";
open (TMP, "<", $tmpFeatureFile);
my %hash;
my %intensityHash;
<TMP>;
while (<TMP>) {
	chomp $_;
	my @data = split(/\t/, $_);
	my ($mz, $ms1ScanNum, $intensity) = ($data[1], $data[2], $data[3]);
	$hash{$ms1ScanNum}{$intensity}{$mz} = $_;
	$intensityHash{$mz} = $intensity;
}
close(TMP);

## Investigate features from one with the highest intensity
## If there's any feature produced from an isotopic peak of another feature,
## those features will be merged 
my %charge;
my %isotope;
foreach my $scan (keys %hash) {
	foreach my $intensity (sort {$b <=> $a} keys %{$hash{$scan}}) {
		foreach my $mz (keys %{$hash{$scan}{$intensity}}) {	
			my $min = 0;
			($scan - 50) < 0 ? $min = 0 : $min = ($scan - 50);				
			for (my $i = $min; $i < ($scan + 50); $i++) {
				my $chargeHash = findCharge($mz, \%{$hash{$i}});
				foreach my $chargedMz (keys %$chargeHash) {
					$charge{$scan}{$mz} = $$chargeHash{$chargedMz};
					if ($intensityHash{$mz} > $intensityHash{$chargedMz}) {
						$isotope{$i}{$chargedMz} = $$chargeHash{$chargedMz};				
					} else {
						$isotope{$scan}{$mz} = $$chargeHash{$mz};
					}
				}
			}
		}
	}
}

open (XML, "<", $mzXML) || die "Cannot open .mzXML file";
my $indexOffset = getIndexOffset(*XML);
my ($indexArray, $lastScan) = getIndexArray(*XML, $indexOffset);
my %scan2RT;
my %RT2scan;
foreach my $index (@$indexArray) {
	my ($scan) = getScanNum(*XML, $index);	
	my ($rt) = getRT(*XML, $index);
	$scan2RT{$scan} = $rt;
	$RT2scan{$rt} = $scan;
}
close (XML);

my $featureFile = $tmpFeatureFile;
$featureFile =~ s/\.tmp\.feature/\.feature/;
open (FEATURE, ">", $featureFile);
print FEATURE "index", "\t", "m\/z", "\t", "z", "\t", "MS1 scan#", "\t", "RT", "\t", 
				"min RT", "\t", "max RT", "\t", "Intensity", "\t", "S\/N", "\t", "Percentage of TF", "\n";
my $index = 0;
foreach my $scan (sort {$a <=> $b} keys %hash) {
	foreach my $intensity (sort {$b <=> $a} keys %{$hash{$scan}}) {
		foreach my $mz (sort {$a <=> $b} keys %{$hash{$scan}{$intensity}}) {	
			next if (defined($isotope{$scan}{$mz}));
			my @data = split("\t", $hash{$scan}{$intensity}{$mz});
			my ($SN, $minRT, $maxRT, $pctTruePeaks) = ($data[4], $data[6], $data[7], $data[8]);
			$index++;
			if (defined($charge{$scan}{$mz})) {
				print FEATURE $index, "\t", sprintf("%.12f", $mz), "\t", $charge{$scan}{$mz}, "\t", $scan, "\t",
								sprintf("%.1f", $scan2RT{$scan}), "\t", $RT2scan{$minRT}, "\t", $RT2scan{$maxRT}, "\t",
								sprintf("%.0f", $intensity), "\t", sprintf("%.1f", $SN), "\t", sprintf("%.4f", $pctTruePeaks), "\n";
			} else {
				print FEATURE $index, "\t", sprintf("%.12f", $mz), "\t", "0", "\t", $scan, "\t",
								sprintf("%.1f", $scan2RT{$scan}), "\t", $RT2scan{$minRT}, "\t", $RT2scan{$maxRT}, "\t",
								sprintf("%.0f", $intensity), "\t", sprintf("%.1f", $SN), "\t", sprintf("%.4f", $pctTruePeaks), "\n";
			}
		}
	}
}
close (FEATURE);

sub findCharge {
	my ($selectMz, $hash) = @_;
	my $C = 1.00335;
	my $intraPpm = 10;	## Decharge ppm
	my $maxCharge = 6;
	my ($lL, $uL) =  ($selectMz, $selectMz + $C + $selectMz * $intraPpm / 1e6); 
	my %chargeHash;	
	foreach my $intensity (sort {$b <=> $a} keys %$hash) {
		foreach my $mz (keys %{$$hash{$intensity}}) {
		# search the previous peak (only one peak) 
			if ($mz > $lL && $mz < $uL) {
				my $diff = 1 / abs($mz - $selectMz);
				my $roundDiff = sprintf("%.0f", $diff);
				next if ($roundDiff == 0 || $roundDiff > $maxCharge);				
				my $var = abs(abs($mz - $selectMz) - ($C / $roundDiff));
				next if ($var > ($intraPpm / 1e6) * $selectMz);   				
				my $charge = $roundDiff;
				$chargeHash{$mz} = $charge;
				return \%chargeHash;
			}
		}
	}
	return \%chargeHash;
}

sub getIndexOffset{
	(*XML) = @_;
	seek (XML, -120, 2);
	my ($indexOffset);
	while (<XML>) {
		next if (!/<indexOffset>/);
		chomp;
		$indexOffset = $_;
		$indexOffset =~ s/\s+<indexOffset>(\d+)<\/indexOffset>.*/$1/o;
		last;
	}
	return $indexOffset;
}

sub getIndexArray {	
	(*XML, my $indexOffset) = @_;
	my @indexArray;
	my $lastScan = 0;
	seek (XML, $indexOffset, 0);
	while (<XML>) {
		next if (/^\s+<index/);
		last if (/^\s+<\/index>.*/);
		chomp;
		next if (/scan/);
		$lastScan++;
		my $index = $_;
		$index =~ s/[\s\t]+\<offset id="\d+"[\s]*>(\d+)<\/offset>.*/$1/o;
		push (@indexArray, $index);
	}
	return (\@indexArray, $lastScan);
}

sub getRT {
	(*XML, my $scanIndex) = @_;
	my $rt;
	seek (XML, $scanIndex, 0);
	while (<XML>) {
		next if (!/retentionTime=/);
		chomp;
		$rt = $_;
		$rt =~ s/.*retentionTime="PT([\d+\.]+)S".*/$1/o;
		last;
	}
	return $rt;
}

sub getScanNum {
	(*XML, my $scanIndex) = @_;
	seek (XML, $scanIndex, 0);
	my $scanNum;
	while (<XML>){
		next if (!/<scan\snum=\"\d+\"/);
		$_ =~ s/^M//g;
		chomp;
		$scanNum = $_;
		$scanNum =~ s/<scan\snum=\"(\d+)\"//;
		$scanNum = $1;
		last;
	}
	return $scanNum;
}
