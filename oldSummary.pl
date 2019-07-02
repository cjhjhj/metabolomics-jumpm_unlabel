#!/usr/bin/perl

use strict;
use warnings;
use lib "/data1/pipeline/dev/JUMPm_v1.8.1.1/";
use Spiders::Params;
use Spiders::Summary;
use Spiders::Job;
use Spiders::Path;

my $library = "/data1/pipeline/dev/JUMPm_v1.8.1.1/";
my $parameter = "jumpm_v1.8.1.1_neg.params";
my $p = Spiders::Params->new('-path' => $parameter);
my $params = $p->parse_param();
my $log_temp = time();
my $LOG;
open ($LOG, ">>", "./$log_temp");
my $dta_path = "/home/htan/2019/IROA_library/neg/align_test/align_test.1/";

## Summarize DB search results
my $job = new Spiders::Job;
$job->set_dta_path($dta_path);
$job->set_parameter($params);	
$job->set_library_path($library);
$job->set_log_file($LOG);	
my $summary = new Spiders::Summary();	
$summary->set_log_file($LOG);
$summary->set_parameter($params);	
my $ms_hash_mol = $summary->missile_summary_unlabel($dta_path);
my ($file_array, $MS1_MS2_matched, $ms2_scan_num) = $job->launch_ms2_jobs_unlabel($ms_hash_mol, $dta_path, $parameter);

#####################
## Finally summary ##
#####################
my $missile_structure = final_summary($ms_hash_mol, $MS1_MS2_matched, $dta_path);


sub final_summary {
	my ($ms_hash_mol, $MS1_MS2_matched, $dta_path) = @_;
	my $missile_structure_with_score = 0;
	my $missile_structure_without_score = 0;
	my $MS2_structure = $dta_path . ".MS2" . ".structure";

	open(OUTPUT,">$dta_path.spectrum_matches") || die "can not open the summary file"; 	
	open(OUTPUTBEST,">$MS2_structure") || die "can not open the summary file";
	print OUTPUT "Index\tMS1 scan\tC12 mass\tN15 mass\tC13 mass\tC12 intensity\tN15 intensity\tC13 intensity\tPscore\tFormula\tAdduct\tTarget-Decoy\tMS2 scan\tName\tSMILES\tInChIKey\tMscore\tPrecursor Type\n";
	print OUTPUTBEST "MS1 scan\tC12 mass\tN15 mass\tC13 mass\tC12 intensity\tN15 intensity\tC13 intensity\tPscore\tFormula\tAdduct\tTarget-Decoy\tMS2 scan\tName\tSMILES\tInChIKey\tMscore\tPrecursor Type\n";

	my $index = 0;
	my $best_index = 0;
	foreach my $scan_missile (sort {$a<=>$b} keys %$ms_hash_mol) {
		foreach my $mz_missile (keys %{$ms_hash_mol->{$scan_missile}}) {
			my $besthit = 1;
			next if ($mz_missile == 0);
			my $flag = 0;
			my %score_hash;		
			foreach my $MS2_scan (keys %{$MS1_MS2_matched->{$scan_missile}->{$mz_missile}}) {
				my $score_file = $MS2_scan . ".MS2.score";
				print "  summarizing scoring file: $score_file \r";
				next if(!-e("$dta_path/$score_file"));
				$flag = 1;						
				open(SCORE,"$dta_path/$score_file") || die "  can not open the file: $score_file \n";
				while(<SCORE>) {
					chomp $_;			
					my @lines = split(/\t/,$_);
					next if(!defined($lines[1]));					
					$lines[1]=~s/\\//g;
					next if(!defined($lines[4]));	
					$score_hash{$lines[4]}{'MS2_scan'} = $MS2_scan ;
					push (@{$score_hash{$lines[4]}{'smiles'}}, $_);					
				}
				close(SCORE);
			}
				
			foreach my $pvalue (sort {$b<=>$a} keys %score_hash) {
				foreach (@{$score_hash{$pvalue}{'smiles'}}) {
					my $score_smile = $_;
					my @score_line=split(/\t/,$score_smile);
					my $ms2_scan = $score_line[0];			
					my $mscore= $score_line[4];
					my $prec_type = $score_line[5];
					my $rank_score= $score_line[6];
					$score_line[1]=~s/\\\(/\(/g;
					$score_line[1]=~s/\\\)/\)/g;
					my $smout = $ms_hash_mol->{$scan_missile}->{$mz_missile}->{$score_line[1]};				
					my ($smout_value) = values %{$smout};
					next if(!defined($smout_value));
					my @smout_line = split(/\t/,$smout_value);
					if (!defined($smout_line[4])) {
							next;
					}
					my $scan = basename($smout_line[4]);
					my $N15_intensity = 0;
					my $C13_intensity = 0;
					if($smout_line[2]!=0) {
						$N15_intensity =$smout_line[6];
					}
					if(defined($smout_line[7]))	{
						$C13_intensity =$smout_line[7];
					}						
					$scan =~s/(\d+)\.MS1/$1/;
					$scan =~s/\.iso//;	
					$index++;
					print OUTPUT $index,"\t",$scan,"\t",sprintf("%.4f",$smout_line[1]),"\t",sprintf("%.4f",$smout_line[2]),"\t",sprintf("%.4f",$smout_line[3]),"\t",$smout_line[5],"\t",sprintf("%.0f",$N15_intensity),"\t",sprintf("%.0f",$smout_line[7]),"\t",sprintf("%.2f",$smout_line[8]),"\t",$smout_line[9],"\t",$smout_line[10],"\t",$smout_line[14],"\t",basename($score_hash{$pvalue}{'MS2_scan'}),"\t",$smout_line[11],"\t",$smout_line[12],"\t",$smout_line[13],"\t",sprintf("%.2f",$pvalue),"\t",$prec_type,"\n";	

					if($besthit==1)	{
						$missile_structure_with_score++;
						$best_index++;
						print OUTPUTBEST $scan,"\t",sprintf("%.4f",$smout_line[1]),"\t",sprintf("%.4f",$smout_line[2]),"\t",sprintf("%.4f",$smout_line[3]),"\t",$smout_line[5],"\t",sprintf("%.0f",$N15_intensity),"\t",sprintf("%.0f",$smout_line[7]),"\t",sprintf("%.2f",$smout_line[8]),"\t",$smout_line[9],"\t",$smout_line[10],"\t",$smout_line[14],"\t",basename($score_hash{$pvalue}{'MS2_scan'}),"\t",$smout_line[11],"\t",$smout_line[12],"\t",$smout_line[13],"\t",sprintf("%.2f",$pvalue),"\t",$prec_type,"\n";

						$besthit=0;						
					}						
				}
			}				
		}
	}
	
	close(OUTPUT);
	close(OUTPUTBEST);
}