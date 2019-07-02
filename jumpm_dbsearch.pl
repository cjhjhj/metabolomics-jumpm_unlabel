#!/usr/bin/perl  

# Load all modules required by JUMPm 
## FindBin for getting current working directory
use FindBin qw($Bin);
use lib "$Bin";
use Parallel::ForkManager;
use Spiders::Params;
use Spiders::ProcessingRAW;
use Spiders::ProcessingMzXML;
use Spiders::LSFQueue;
use Spiders::SGEQueue;
use Spiders::Path;
use Getopt::Long;
use Cwd;
use Cwd 'abs_path';
use File::Basename;

####################
## Initialization ##
####################
my $library = $Bin;
my $RLibrary = $Bin . "/R/";
if (!-e($RLibrary)) {
	print "R scripts are NOT accessible from jump -m\n";
	exit;
}
my ($help, $paramFile, $rawFile);
GetOptions('-help|h'=> \$help, 
			'-p=s' => \$paramFile,
			);
if (!-e ($paramFile)) {
	print "Please input the parameter file\n\n";
	exit;
}
$paramFile = abs_path($paramFile);
usage() if ($help || !defined($paramFile));

## Define a log file usings time() function
my $tmpLog = time();
my $LOG;
## For unlabeled data, 
## the log file should be open with append option to include R script captions
open ($LOG, ">>", "./$tmpLog");
print "\n\n  Initializing JUMPm program\n\n";
print $LOG "\n\n  Initializing JUMPm program\n\n";
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time); 

## Loading parameters
my $p = Spiders::Params -> new('-path' => $paramFile);
my $params = $p -> parseParams();

my ($ms2Path) = @ARGV; 

#########################################################################
## Database search to identify metabolites using individual .MS2 files ##
#########################################################################
my $queue;
if ($$params{'cluster'} == 1 && ($$params{'cluster_type'} eq "LSF" || $$params{'job_management_system'} eq "LSF")) {
	$queue = Spiders::LSFQueue->new();
} elsif ($$params{'cluster'} == 1 && ($$params{'cluster_type'} eq "SGE" || $$params{'job_management_system'} eq "SGE")) {
    $queue = Spiders::SGEQueue->new();
}
my @ms2FileArray = glob("$ms2Path/*.MS2");
my $nFiles = scalar(@ms2FileArray);
my $maxJobs = 200;
my $filesPerJob;
if ($nFiles <= $maxJobs) {
	$maxJobs = $nFiles;
	$filesPerJob = 1;
} else {
	$filesPerJob = int($nFiles / $maxJobs) + 1;
}
my $nTotalJobs = int($nFiles / $filesPerJob - 0.0001) + 1;
my $nJobs = 0;
my $randNum = int(rand(100));
my %jobIDs;
for (my $i = 0; $i < $nTotalJobs; $i++) {
	$nJobs++;
	my $jobName = "sch_${randNum}_$nJobs";
	my $command = "";
	for (my $j = 0; $j < $filesPerJob; $j++) {
		my $k = $filesPerJob * $i + $j;
		last if ($k >= $nFiles);
		$command .= "perl $Bin/databaseSearch.pl $paramFile $ms2FileArray[$k]\n";
	}
	my $job = $queue -> submit_job($ms2Path, $jobName, $command);
	$jobIDs{$job} = 1;
	print "\r  $nJobs database search jobs are submitted";
}
print "\n  You submitted $nJobs jobs for database search\n";
checkJobStatus($nJobs, \%jobIDs, $queue);

sub checkJobStatus {
	my ($nJobs, $jobIDs, $queue) = @_;
	my $nFinishedJobs = 0;
	my $jobInfo = 1;
	my ($username) = getpwuid($<);
	$| = 1;
	while($jobInfo) {
		my @jobStatusArray = @{$queue->get_running_jobs($jobIDs)};
		if (@jobStatusArray) {  # i.e. There are some jobs in the queue (may include other jobs like searching)
			my @jobIDsInQueue;
			foreach (@jobStatusArray) {
				my ($jobID) = ($_ =~ /([0-9]+)\s*/);
				if (defined $$jobIDs{$jobID}) {
					push (@jobIDsInQueue, $jobID);
				}
			}
			if (@jobIDsInQueue) { # i.e. There are jobs of interest in the queue
				$nFinishedJobs = $nJobs - scalar(@jobIDsInQueue);
				print "\r  $nFinishedJobs jobs are finished";
				if ($nFinishedJobs == $nJobs) {
					$jobInfo = 0;
				} else {
					$jobInfo = 1;
				}
			} else {        # i.e. There are jobs in the queue, but all jobs of interest are finished
				print "\r  $nJobs jobs are finished";
				$jobInfo = 0;
			}
		} else {        # i.e. There's no job in the queue
			print "\r  $nJobs jobs are finished";
			$jobInfo = 0;
		}
		sleep(5);
	}
	$| = 0;
}


## ============================================================================= ##

#######################
## Feature alignment ##
#######################

## Using @featureFileArray


###############################################################
## Generation of (consolidated) MS2 spectrum of each feature ##
###############################################################

## Directory handling should be carefully treated


#############################################
## Database search to identify metabolites ##
#############################################

## Fully-aligned feature table has m/z values of features at the 1st column

## 1. First search; against formula and structure database using m/z only

## 2. Second search; comparison between observed (consolidated) MS2 and theoretical MS2 (produced from fragmentation)




=head

# Get the R library 
## Define the R version as R-3.1.0
my $library = $Bin;
my $R_library = $Bin . "/R-3.1.0/";
if (!-e($R_library)) {
	print "\n\n  Please build the R (version: 3.1.0) under $library \n\n";
	exit;
}

# Define JUMPm version
my $VERSION = "1.8.1";
my $progname = $0;
$progname =~ s@(.*)/@@i;

my ($help, $parameter, $raw_file);
GetOptions('-help|h'=>\$help,
			'-p=s'=>\$parameter,
		);
if (!-e ($parameter)) {
	print "please input the paramter file\n\n";
	exit;
}
usage() if ($help || !defined($parameter));

# Define log file as time
my $log_temp = time();
my $LOG;
if ($params->{'labeled_data'} == 1) {
	open($LOG, ">./$log_temp");
} elsif ($params->{'labeled_data'} == 0) {
	## For unlabeled data, the log file should be open with append option to include R script captions
	open ($LOG, ">>", "./$log_temp");
}
print "\n\n  Initializing JUMPm program (version $VERSION)\n\n";
print $LOG "\n\n  Initializing JUMPm program (version $VERSION)\n\n";
my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time); 

my $p = Spiders::Params->new('-path'=>$parameter);
my $params = $p->parse_param();
$p->parameter_check($params);

# searching a spectral library to the database
if ($params->{'labeled_data'} == 2) {
	my $dta_path = $ARGV[0];
	my @file_array = glob("$dta_path/*.dta");
	my $job = new Spiders::Job;	
	$job->set_dta_path($dta_path);
	$job->set_parameter($params);	
	$job->set_library_path($library);
	$job->set_log_file($LOG);	
	$job->create_mscript_unlabel_dta();
	$random = int(rand(100));
	$parameter = abs_path($parameter);			
	$job->runjobs(\@file_array,$dta_path, "sch_${random}", 0, "Undef", $parameter);	
	my $summary = new Spiders::Summary();	
	$summary->set_log_file($LOG);	
	my $ms_hash_mol = $summary->missile_summary_unlabel($dta_path);
	my @frag_job_list=();
	my $prec_type="C12";
	my $target_decoy = "target";
	my $MS1_MS2_matched;
	foreach my $dta (@file_array) {
		my ($filename, $dirs, $suffix) = fileparse($dta);
		foreach $mass (keys %{$ms_hash_mol->{$filename}}) {
			$MS1_MS2_matched->{$filename}->{$mass}->{$filename} = 1;
			foreach $smile (keys %{$ms_hash_mol->{$filename}->{$mass}}) {
				$smile =~ s/\(/\\\(/g;
				$smile =~ s/\)/\\\)/g;
				push (@frag_job_list, "perl $library/frag_shell.pl -dtafile $dta  -smile \"$smile\" -mass 50 -depth 2 -param  $parameter -ptype $prec_type -target_decoy $target_decoy");
			}
		}
	}
	$random = int(rand(100));
	$job->runjobs(\@frag_job_list, "$dta_path", "frag_${random}", 0, 0, $parameter);
	my $missile_structure = $summary->final_summary_dta($ms_hash_mol, $MS1_MS2_matched, $dta_path);
	exit;
}

## Create the path for multiple raw files
my %rawfile_hash;
print "  Using the following rawfiles:\n";
print $LOG "  Using the following rawfiles:\n";
foreach $arg (sort @ARGV) {
	my @suffixlist = ();
	push @suffixlist, ".raw";
	push @suffixlist, ".RAW";
	push @suffixlist, ".mzXML";

	if ($arg =~ /.[raw|RAW|mzXML]/) {
		print "  $arg","\n";
		print $LOG "  $arg","\n";
	}
}

## Start to process each dataset
## prepare an array of raw data
my @fileArray; 
foreach $arg (sort @ARGV) {
	my @suffixlist = ();
	push @suffixlist,".raw";
	push @suffixlist,".RAW";
	push @suffixlist,".mzXML";

	if ($arg =~ /.[raw|RAW|mzXML]/) {
		push (@fileArray, $arg);
		my ($filename, $directory, $suffix) = fileparse($arg, @suffixlist);
		system(qq(mkdir $directory/$filename >/dev/null 2>&1));
		system(qq(mv $arg $directory/$filename >/dev/null 2>&1));
		my $datafile = "$directory/$filename";
		if (!-e($datafile)) {
			print "  There is no such raw file ($arg)!\n";
			exit;
		}		
		my $path = new Spiders::Path($datafile);
		my $list = $path->make_directory_list();
		if (@$list) {
			$newdir = $path->choose_dir_entry($list, "  Choose a .out file directory", $newdir);
		} else {
			$newdir	= $filename . ".1";
		}
		print "  Using: $newdir\n";
		print $LOG "  Using: $newdir\n";		
		$path->add_subdir($newdir);
		my $dir =  $path->basedir() . "/$newdir";
		my $rawfile = "$datafile/$arg";
		$rawfile_hash{$rawfile} = $dir; 
	}
}

##############################
## Jump -m for LABELED data ##
##############################
if ($params->{'labeled_data'} == 1) {
## start processing each dataset
	foreach my $raw_file (sort keys %rawfile_hash) {		
		print "\n  Searching data: $raw_file\n";
		print $LOG "\n  Searching data: $raw_file\n";	

		## Get system time 
		my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time); 
		print "  Start: ";
		print $LOG "  Start: ";	
		printf "%4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;	
		printf $LOG "%4d-%02d-%02d %02d:%02d:%02d\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;	
		
		my $curr_dir = getcwd;
		my $dta_path = $rawfile_hash{$raw_file};
		if ($raw_file =~/\.mzXML/) {
			$raw_file =~ s/\.mzXML//g;
		}

		## convert raw data into mzXML
		my $proc_raw = new Spiders::ProcessingRAW();
		$proc_raw->set_log_file($LOG);		
		$proc_raw->set_raw_file($raw_file);
		print "  Converting .raw into .mzXML file\n";
		print $LOG "  Converting .raw into .mzXML file\n";
		my $proc_xml = new Spiders::ProcessingMzXML();	
		$proc_xml ->set_dta_path($dta_path);
		$proc_xml ->set_log_file($LOG);	
		my (%ms_hash, %msms_hash, @mz_array);		
		my $MS1_peaks_sequenced;
		my ($msms_hash_corrected, $mz_array_corrected, $meanMassShift, $stdMassShift) = (0, 0, 0, 0);

		## if the data is centroid mode 	
		if ($params->{data_acquisition_mode} == 1) {	
			my $mzXML = $proc_raw->raw2mzXML("centroid");
			print "\n\n  Using centroid mode\n";
			print $LOG "\n\n  Using centroid mode\n";		
			print "  Extracting MS1 and MS2 peaks from .mzXML\n";
			print $LOG "  Extracting MS1 and MS2 peaks from .mzXML\n";	
			$proc_xml ->set_mzXML_file($mzXML);

			## preprocessing
			$proc_xml->set_parameter($params);
			$proc_xml->generate_hash_dta(\%ms_hash, \%msms_hash, \@mz_array, $params);
			my $ms1N = scalar(keys %{$ms_hash{'surveyhash'}});
			my $ms2N = scalar(keys %msms_hash)-scalar(keys %{$ms_hash{'surveyhash'}});
			printf("\n  There are %d MS and %d MS/MS in the entire run\n", $ms1N , $ms2N);
			printf $LOG ("\n  There are %d MS and %d MS/MS in the entire run\n", $ms1N , $ms2N);

			if ($params->{mass_correction}) {
				print "\n\n  Mass correction\n";
				print $LOG "\n\n  Mass correction\n";	
				my $masscorr = new Spiders::MassCorrection();
				$masscorr->set_log_file($LOG);
				## No correction for MS1, but it will perform correction at the database search step
				($msms_hash_corrected, $mz_array_corrected, $meanMassShift, $stdMassShift) = $masscorr->massCorrection(\%ms_hash, \%msms_hash, \@mz_array, $params);
				printf $LOG ("  Calculated mass-shift: mean = %.5f ppm and SD = %.5f ppm\n", $meanMassShift, $stdMassShift); 
				%msms_hash = %$msms_hash_corrected;
				@mz_array = @$mz_array_corrected;	
			}
			
			my $decharge = new Spiders::Decharge();
			$decharge->set_parameter($params);
			$decharge->decharge(\%ms_hash, \%msms_hash, \@mz_array);

			print "  Generating MS1 dta files\n";	
			print $LOG "  Generating MS1 dta files\n";			
			$proc_xml->generate_MS1_file(\@mz_array);
			if($ms2N>0) {
				$MS1_peaks_sequenced = $proc_xml->generate_MS2_file(\%msms_hash, \@mz_array);
			}
		}

		## if the data is profile mode 
		elsif($params->{data_acquisition_mode} == 2) {
			my $mzXML = $proc_raw->raw2mzXML("profile");
			print "\n\n  Using profile mode\n";	
			print $LOG "\n\n  Using profile mode\n";		
			print "  Feature detection\n";
			print $LOG "  Feature detection\n";	

			## Extract MS1 and MS2 data
			$proc_xml ->set_mzXML_file($mzXML);
			$proc_xml->set_parameter($params);
			print $LOG "  Gathering scan information\n";			
			$proc_xml->generate_hash_dta(\%ms_hash, \%msms_hash, \@mz_array, $params);
			my $ms_hash_peak_detection = $proc_xml->{'_mshash4peakdetection'};
			my $ms1N = scalar(keys %{$ms_hash{'surveyhash'}});
			my $ms2N = scalar(keys %msms_hash)-scalar(keys %{$ms_hash{'surveyhash'}});
			
			printf("\n  There are %d MS and %d MS/MS in the entire run\n", $ms1N , $ms2N);
			printf $LOG ("\n  There are %d MS and %d MS/MS in the entire run\n", $ms1N , $ms2N);

			## Feature detection for MS1 peaks ######################################
			my $peakdetection = new Spiders::PeakDetection();
			$peakdetection->set_log_file($LOG);
			my ($ms, $peak) =$peakdetection->Gen_3Dpeaks($params, $mzXML, $dta_path);
			## Feature detection end ################################################

			my $correction_raw = $peakdetection->generate_correction_file($peak, $dta_path);	

			## if it is a positive mode data, and would like to perform mass correction	
			## we can only perform correction for positive data using 445 peak
			## we can not detect 443 peak in negative data ???
	 		if ($params->{mass_correction} and $params->{mode} == 1) {
				## Mass correction start #############################################
				print "\n  Mass correction\n";
				print $LOG "\n  Mass correction\n";	
				my $masscorr = new Spiders::MassCorrection();
				my $mass_shift_file = $dta_path . ".mass_shift";
				$masscorr->set_log_file($LOG);
				($ms_hash_corrected, $msms_hash_corrected, $mz_array_corrected, $mass_shift_mean_all) = $masscorr->massCorrection(\%ms_hash, \%msms_hash, \@mz_array,$correction_raw, $params, $mass_shift_file);
				## save all mass shift under a hiden folder .mass_shift
				store $mass_shift_mean_all, "$dta_path/.mass_shift";
				@mz_array = @$mz_array_corrected;
				## Mass correction end ################################################

				## Generate MS1 and MS2 peak files
				#$peakdetection->generate_MS1_file($ms, $peak, $dta_path);
				
				## Taking action: Mass correction for MS1 data
				$peakdetection->generate_MS1_file_based_feature($peak, $dta_path, $params);
				$MS1_peaks_sequenced = $proc_xml->generate_MS2_file($msms_hash_corrected, \@mz_array);
			} else {
				#$peakdetection->generate_MS1_file($ms, $peak, $dta_path);
				$peakdetection->generate_MS1_file_based_feature($peak, $dta_path, $params);
				if ($ms2N > 0) {			
					$MS1_peaks_sequenced = $proc_xml->generate_MS2_file(\%msms_hash, \@mz_array);
				}
			}
			print "\n";
			print $LOG "\n";		
		} else {
			system(qq(cp $log_temp $dta_path/JUMPm.log));		
			print " please set the right data_acquisition_mode\n";
			exit;
		}
			
		my @dtafiles = glob("$dta_path/*.MS1");
		my $dta = new Spiders::Dta; 	
		my $job = new Spiders::Job;
		$job->set_log_file($LOG);
		$job->set_library_path($library);
		$job->set_R_library_path($R_library);		
		$job->set_dta_path($dta_path);
		$job->set_parameter($params);
		
		my $summary = new Spiders::Summary();	
		$summary->set_library_path($R_library);
		$summary->set_log_file($LOG);
		$summary->set_parameter($params);	
		
		if ($params->{data_acquisition_mode} == 2) {	
			$summary->generate_feature_table($dta_path);
		}
		
		print "  You are searching a labeled data set\n";
		print $LOG "  You are searching a labeled data set\n";
		if ($params->{'labeled_data'} == 1) {		
			$job->create_pscript;
		} elsif ($params->{'labeled_data'} == 2) {		
			$job->create_mscript_partial;
		}
		my $random = int(rand(100));
		my $ratio;
		my $defect;
		$ratio->{'NC'} = 0;
		$ratio->{'CC'} = 0;
		$parameter = abs_path($parameter);			
		if ($params->{loading_normalization}) {		
			print "  Estimating the labeling mix ratio (it takes time and please be patient!)\n";
			print $LOG "  Estimating the labeling mix ratio (it takes time and please be patient!)\n";
			$job->runjobs(\@dtafiles, $dta_path, "pair_${random}", 0, 0, $parameter);		
			($ratio,$defect) = $summary->calculate_intensity_ratio($dta_path);
			my $NC_ratio = 2**$ratio->{'NC'};
			my $CC_ratio = 2**$ratio->{'CC'};
			print "  Intensity ratio between N15 and C12 is $NC_ratio \n";
			print "  Intensity ratio between C13 and C12 is $CC_ratio \n";
			print $LOG "  Intensity ratio between N15 and C12 is $NC_ratio \n";
			print $LOG "  Intensity ratio between C13 and C12 is $CC_ratio  \n";
			print "  Mass defect between N15 and C12 is $defect->{'NC'}->[0] \n";
			print "  Mass defect between C13 and C12 is $defect->{'CC'}->[0] \n";
			print $LOG "  Mass defect between N15 and C12 is $defect->{'NC'}->[0] \n";
			print $LOG "  Mass defect between C13 and C12 is $defect->{'CC'}->[0] \n";
		} else {
			print "  No loading normalization\n";
			print $LOG "  Not loading normalization\n";
			$ratio->{'CC_std'} = 1;			
			$ratio->{'NC_std'} = 1;			
			$defect->{'NC'}->[0] = 1.00335;
			$defect->{'CC'}->[0] = 0.99703;
			$defect->{'NC'}->[1] = 1;
			$defect->{'CC'}->[1] = 1;			
		}
		print "\n  Pair scoring, searching formulas and structures (it takes time and please be patient!)\n";
		print $LOG "\n  Pair scoring, searching formulas and structures (it takes time and please be patient!)\n";
		if ($params->{'labeled_data'} == 1) {		
			$job->create_mscript;
		} elsif ($params->{'labeled_data'} == 2) {
			system(qq(mkdir $dta_path/MS1));
			system(qq(cp $dta_path/*.MS1 $dta_path/MS1/));
			$job->create_mscript_partial2;
		}
		$random = int(rand(100));
				
		$job->runjobs(\@dtafiles, $dta_path, "sch_${random}", $ratio, $defect, $parameter);
		print "\n\n  Summarizing pair results\n";
		print $LOG "\n\n  Summarizing pair results\n";
		$summary->pair_summary($dta_path);
		print "\n  Summarizing missile results";
		print $LOG "\n  Summarizing missile results";	
		my $ms_hash_mol = $summary->missile_summary($dta_path);
		my $missile_scan = scalar keys (%$ms_hash_mol);
			
		if ($missile_scan > 0) {	
			print "\n  Performing MS2 matching and scoring\n\n";
			print $LOG "\n  Performing MS2 matching and scoring\n\n";	
			my ($file_array, $MS1_MS2_matched, $ms2_scan_num) = $job->luanch_ms2_jobs(\%msms_hash, $ms_hash_mol, $dta_path, $parameter);
			print "\n  $ms2_scan_num MS2 scans match to MISSILES (one missile may have multiple MS2 matches)\n";
			print $LOG "\n  $ms2_scan_num MS2 scans match to MISSILES (one missile may have multiple MS2 matches)\n";
			if ($ms2_scan_num > 0) {
				$random = int(rand(100));
				$job->runjobs($file_array, $dta_path, "frag_${random}", 0, 0, $parameter);
				
				## Finally summary
				my $missile_structure = $summary->final_summary($ms_hash_mol, $MS1_MS2_matched, $dta_path);
			}
		}

		($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time); 
		print "  JUMPm finished: ";
		print $LOG "  JUMPm finished: ";	
		printf "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;	
		printf $LOG "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;

		system(qq(mkdir $dta_path.intermediate));
		if ($params->{'labeled_data'} == 1) {
			system(qq(mv $dta_path.pair.all.comb $dta_path.intermediate/));	
			system(qq(mv $dta_path.pair.all.comb.uniq $dta_path.intermediate/));
		}
		system(qq(mv $dta_path.tmp.feature $dta_path.intermediate/));
		system(qq(mv $dta_path.missile $dta_path.intermediate/));
		system(qq(mv $dta_path.missile.unique.pass $dta_path.intermediate/));
		system(qq(mv $dta_path.MS2.structure $dta_path.intermediate/));
		system(qq(cp $parameter $dta_path.intermediate/));

		my $dir = dirname($log_temp);
		system(qq(cp -rf $dir/JUMPm.log $dta_path.intermediate/));
		system(qq(mv $log_temp $dir/JUMPm.log));	
	}
	# end of processing LABELED data

################################
## Jump -m for UNLABELED data ##
################################
} elsif ($params->{'labeled_data'} == 0) {
	########################
	## .raw file handling ##
	########################
	## Note: in this script, it is supposed that "data_acquisition_mode" is always 2
	my $random = int(rand(100));
	my @mzxml_files = ();
	my @mzxml_paths = ();
	my $job = new Spiders::Job;
	$job->set_log_file($LOG);
	$job->set_library_path($library);
	$job->set_parameter($params);
	foreach my $raw_file(sort keys %rawfile_hash) {
		my $dta_path = $rawfile_hash{$raw_file};	# $dta_path should be defined before removing ".mzXML" from $raw_file
		if ($raw_file =~ /\.mzXML/) {
			$raw_file =~ s/\.mzXML//g;
		}
		## Convert raw data into mzXML
		my $proc_raw = new Spiders::ProcessingRAW();
		$proc_raw->set_log_file($LOG);		
		$proc_raw->set_raw_file($raw_file);
		$mzXML = $proc_raw->raw2mzXML("profile");

		## Job for feature detection
		my $mzxml_path = dirname(abs_path($mzXML));
		$job->set_mzxml_path($mzxml_path);
		$job->set_dta_path($dta_path);
		$job->create_fscript_unlabel;
		push (@mzxml_files, basename($mzXML));
		push (@mzxml_paths, $mzxml_path);
		
	}
	print "\n  Feature detection is being performed\n";
	print "  Please be patient. It may take a couple of hours\n";
	$job->runjobs(\@mzxml_files, \@mzxml_paths, "feat_${random}", 0, 0, abs_path($parameter));

	#######################
	## Feature alignment ##
	#######################
	## Create alignment.params file
	open (AP, ">", "alignment.params");
	foreach my $key (sort {$a cmp $b} keys %$params) {
		if ($key eq "tol_initial") {
			print AP "$key" . " = " . $$params{$key} . "\n";
		} elsif ($key eq "sd_width") {
			print AP "$key" . " = " . $$params{$key} . "\n";
		} elsif ($key eq "rescue") {
			print AP "$key" . " = " . $$params{$key} . "\n";
		#} elsif ($key eq "sd_width_rescue") {
		#	print AP "$key" . " = " . $$params{$key} . "\n";
		#} elsif ($key eq "intensity_level_rescue") {
		#	print AP "$key" . " = " . $$params{$key} . "\n";
		#} elsif ($key eq "pct_intensity_ratio_rescue") {
		#	print AP "$key" . " = " . $$params{$key} . "\n";
		} elsif ($key eq "pct_full_alignment") {
			print AP "$key" . " = " . $$params{$key} . "\n";
		} elsif ($key eq "output_name") {
			print AP "$key" . " = " . $$params{$key} . "\n";
		} elsif ($key eq "rt_tolerance_unit") {
			print AP "$key" . " = " . $$params{$key} . "\n";
		} elsif ($key eq "rt_tolerance_value") {
			print AP "$key" . " = " . $$params{$key} . "\n";
		} elsif ($key eq "mz_tolerance_unit") {
			print AP "$key" . " = " . $$params{$key} . "\n";
		} elsif ($key eq "mz_tolerance_value") {
			print AP "$key" . " = " . $$params{$key} . "\n";
		}
	}
	my @featureFileArray;
	foreach my $key (sort {$a cmp $b} keys %rawfile_hash) {
		my $featureFileName = $rawfile_hash{$key} . ".feature";
		push (@featureFileArray, $featureFileName);
	}
	my $featureFiles = join(",", @featureFileArray);

	## Execute R script for feature alignment
	my $alignDir = "./align_" . $$params{'output_name'};
	system(qq(mkdir $alignDir >/dev/null 2>&1));
	my $command = "Rscript $Bin/R/alignment.R alignment.params $featureFiles $log_temp $Bin $alignDir";
	system ($command);

	#########################################
	## Processing MS2 spectra for features ##
	#########################################
	## Create featureToMs2.params file
	open (FP, ">", "featureToMs2.params");
	my $featureFile = $alignDir . "/" . $$params{'output_name'} . "_fully_aligned.feature";
	print FP "feature_file = $featureFile\n";
	## "tol_isolation" in featureToMs2 should be the same as "isolation_window" parameter
	my $tol_isolation = $$params{'isolation_window'} / 2;
	print FP "tol_isolation = $tol_isolation\n";
	foreach my $key (sort {$a cmp $b} keys %$params) {
	        if ($key eq "tol_precursor") {
	                print FP "$key" . " = " . $$params{$key} . "\n";
	        } elsif ($key eq "tol_intra_ms2_consolidation") {
	                print FP "$key" . " = " . $$params{$key} . "\n";
	        } elsif ($key eq "tol_inter_ms2_consolidation") {
	                print FP "$key" . " = " . $$params{$key} . "\n";
	        }
	}
	for (my $i = 0; $i < scalar(@fileArray); $i++) {
		foreach my $key (keys %rawfile_hash) {
			if ($fileArray[$i] eq basename($key)) {
				my $mzXML = $key;
				$mzXML =~ s/\.(raw|RAW|Raw)/\.mzXML/;
				if ($i == 0) {
					print FP "mzXML = " . $mzXML . "\n";
				} else {
					print FP "$mzXML\n";
				}
			}
		}
	}
	close (FP);

	## Execute R script for processing MS2 spectra
	my ($dta_path) = $newdir =~ /(\.\d+$)/;
	$dta_path = abs_path($alignDir) . "/align_" . $$params{'output_name'} . $dta_path;
	$command = "Rscript $Bin/R/featureToMs2.R featureToMs2.params $dta_path $log_temp";
	system ($command);

	#####################
	## Database search ##
	#####################
	print "\n  Searching formulas and structures (it takes time and please be patient!)\n";		
	print $LOG "\n  Searching formulas and structures (it takes time and please be patient!)\n";
	$job->set_dta_path($dta_path);
	my @dtafiles = glob("$dta_path/*.MS2");
	if ($params->{'library_search'} == 1) {
		$job->create_mscript_unlabel_library();
	} else {
		$job->create_mscript_unlabel();
	}
	$random = int(rand(100));
	$parameter = abs_path($parameter);
	$job->runjobs(\@dtafiles, $dta_path, "sch_${random}", 0, 0, $parameter);

	## Summarize DB search results
	my $summary = new Spiders::Summary();	
	$summary->set_log_file($LOG);
	$summary->set_parameter($params);	
	my $ms_hash_mol = $summary->missile_summary_unlabel($dta_path);
	my ($file_array, $MS1_MS2_matched, $ms2_scan_num) = $job->launch_ms2_jobs_unlabel($ms_hash_mol, $dta_path, $parameter);
	print "\n  $ms2_scan_num MS2 scans match to MISSILES (one missile may have multiple MS2 matches)\n";
	print $LOG "\n  $ms2_scan_num MS2 scans match to MISSILES (one missile may have multiple MS2 matches)\n";
	$random = int(rand(100));
	$job->runjobs($file_array, $dta_path, "frag_${random}", 0, 0, $parameter);

	#####################
	## Finally summary ##
	#####################
	my $missile_structure = $summary->final_summary_unlabel($ms_hash_mol, $MS1_MS2_matched, $dta_path);
	($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time); 
	print "  JUMPm finished: ";
	print $LOG "  JUMPm finished: ";	
	printf "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;	
	printf $LOG "%4d-%02d-%02d %02d:%02d:%02d\n\n", $year + 1900, $mon + 1, $mday, $hour, $min, $sec;
	system(qq(mkdir $dta_path.intermediate));
	system(qq(mv $dta_path.tmp.feature $dta_path.intermediate/));
	system(qq(mv $dta_path.missile $dta_path.intermediate/));
	system(qq(mv $dta_path.missile.unique.pass $dta_path.intermediate/));
	system(qq(mv $dta_path.MS2.structure $dta_path.intermediate/));
	system(qq(cp $parameter $dta_path.intermediate/));
	my $dir = dirname($log_temp);
	system(qq(cp -rf $dir/JUMPm.log $dta_path.intermediate/));
	system(qq(mv $log_temp $dir/JUMPm.log));
}

close($LOG);
	
sub usage {
print <<"EOF";
	
################################################################
#                                                              #
#       **************************************************     # 
#       ****                                          ****     # 
#       ****                 JUMPm                    ****     # 
#       ****        version $VERSION                  ****     # 
#       ****        Xusheng Wang / Junmin Peng        ****     # 
#       ****         Copyright (C) 2012 - 2013        ****     # 
#       ****            All rights reserved           ****     # 
#       ****                                          ****     # 
#       **************************************************     # 
#                                                              #
################################################################

Usage: $progname -p parameterfile rawfile.raw 
	or
       $progname -p parameterfile rawfile.mzXML
	
EOF
exit 1;
}
