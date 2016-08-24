#!/usr/bin/perl

#####################################################################
#Read target regions into hash structure
#Go through each line of pileup file and identify whether it lies within target region - if yes then increment corresponding target counts in hash
#Calculate a range of coverage data
######################################################################

use strict;
use warnings;
use Getopt::Long;

#CALL: Coverage_from_pileup_Jan11.pl --targets "targetsfile.txt"
my $path1="rootfilepath";
my $path2="outfilepath";
my $file="infilenameANDpath";
my $batchID="overall_sample_ID";
my $TargetFile="/users/nhrg/lustre/hg19/CCDS_21Feb2012_exonslist_nonOV.txt_1-basedpos";

my $Results=GetOptions("inPath1=s"=>\$path1, "inPath2=s"=>\$path2, "inFile=s"=>\$file, "targets=s"=>\$TargetFile, "batchID=s"=>\$batchID);

my %files;
my @Patient_ID;
if($file =~ /(\S+\/)(\S+)\.pileup$/){$files{$2}=$file; push(@Patient_ID,$2);}

#open target regions file
open (FILE, $TargetFile) || die "File not found\n";
print "$TargetFile\n";
#Create hash of target regions and hash of no. of entries per chr
my %Targetregions;
my %Targets_per_chr;
my $targetbase_count=0;
while (<FILE>) {
	chomp($_);
	my @targetlinesplit = split(/\t/, $_);
	my $chr=$targetlinesplit[0];
	my $T_start=$targetlinesplit[1];
	my $T_end=$targetlinesplit[2];
	my $T_info=$targetlinesplit[3];
	my $Tlength =($T_end+1)-$T_start;
	$targetbase_count+=$Tlength;
	if(!exists $Targets_per_chr{$chr}){$Targets_per_chr{$chr}=0;}
	$Targetregions{$chr}[$Targets_per_chr{$chr}][0]=$T_start;
	$Targetregions{$chr}[$Targets_per_chr{$chr}][1]=$T_end;
	$Targetregions{$chr}[$Targets_per_chr{$chr}][2]=0; #num. bases covered 1+
	$Targetregions{$chr}[$Targets_per_chr{$chr}][3]=0; #total coverage per region (1+)
	$Targetregions{$chr}[$Targets_per_chr{$chr}][4]=$T_info; #target region name
	$Targetregions{$chr}[$Targets_per_chr{$chr}][5]=0; #num. bases covered 20+
	$Targetregions{$chr}[$Targets_per_chr{$chr}][6]=0; #num. bases covered 5+
	$Targetregions{$chr}[$Targets_per_chr{$chr}][7]=0; #num. bases covered 10+
	$Targetregions{$chr}[$Targets_per_chr{$chr}][8]=0; #num. bases covered 30+
	$Targets_per_chr{$chr}++;
}
close FILE; 

my $Targetregions_ref = \%Targetregions;
my $Targets_per_chr_ref = \%Targets_per_chr;

#open output file for targets covered results
my $outfile = $path1."/Coverage_from_pileup-targets_hg19_".$batchID.".txt";
open(OUT, ">>$outfile") || die "Cannot open file \"$outfile\" to write to!\n";
     
#open output file for bases covered results
my $outfile2 = $path1."/Coverage_from_pileup-bases_on_target_hg19_".$batchID.".txt";
open(OUT2, ">>$outfile2") || die "Cannot open file \"$outfile2\" to write to!\n";

#open output file for mean coverage results
my $outfile4 = $path1."/Coverage_from_pileup-mean_targetbase_coverage_hg19_".$batchID.".txt";
open(OUT4, ">>$outfile4") || die "Cannot open file \"$outfile4\" to write to!\n";
 
#read all pileup files
my $num_patients = scalar @Patient_ID;
my $Results_bases;
my $Results_targets;
my $Results_mean_cov;
for(my $p1=0; $p1<$num_patients; $p1++){
	 ($Results_bases, $Results_targets, $Results_mean_cov) = Read_pileup($path2, $p1, $Patient_ID[$p1],$files{$Patient_ID[$p1]},$Targetregions_ref,$Targets_per_chr_ref, $targetbase_count);
	 print OUT "$Results_targets\n"; 
	 print OUT2 "$Results_bases\n"; 
	 print OUT4 "$Results_mean_cov\n";
	 }
close OUT;
close OUT2;
close OUT4;

##SUBROUTINES###########################################################################
sub Read_pileup{
	my ($path, $p1, $ID, $pileup_file, $ref_Tregions, $ref_TperChr, $Tbase_count)=@_;
	my %T_regions = %$ref_Tregions;
	my %T_perChr = %$ref_TperChr;
	my %Startfrom;
	my $Bcount_30fold=0;
	my $Bcount_20fold=0;
	my $Bcount_10fold=0;
	my $Bcount_5fold=0;
	my $Bcount_1fold=0;
	my $totalcov=0;
	my $max=0;
	my $min=10000000;
	
	open PILEUP, $pileup_file or die "Pileup file not found\n";
	while (<PILEUP>){
		my $line = $_;
		chomp $line;
		my @linesplit = split (/\t/,$line);
		my $chr = $linesplit[0];
		if(!exists $T_regions{$chr}){next;}
		my $pos = $linesplit[1];
		my $cov = $linesplit[3];
 
	#loop through chr specific targets
	if(!exists $Startfrom{$chr}){$Startfrom{$chr}=0;}
	for (my $count = $Startfrom{$chr}; $count<$T_perChr{$chr}; $count++){	
		
		#if pos is less than tstart, make startfrom{chr} equal to count and move onto next line of file
		if($pos < $T_regions{$chr}[$count][0]){$Startfrom{$chr} = $count; last;}
		
		#if pos is greater than or equal to tstart and pos is less than or equal to tend, add cov to hash and make startfrom{chr} equal to count and move onto next file line
		if(($pos >= $T_regions{$chr}[$count][0]) && ($pos <= $T_regions{$chr}[$count][1])){
			if($cov>=30){$Bcount_30fold++; $T_regions{$chr}[$count][8]++;} #count 30fold covered bases on target - total and per target
			if($cov>=20){$Bcount_20fold++; $T_regions{$chr}[$count][5]++;} #count 20fold covered bases on target - total and per target
			if($cov>=10){$Bcount_10fold++; $T_regions{$chr}[$count][7]++;} #count 10fold covered bases on target - total and per target
			if($cov>=5){$Bcount_5fold++; $T_regions{$chr}[$count][6]++;} #count 5fold covered bases on target - total and per target
			$Bcount_1fold++; #count bases on target
			$T_regions{$chr}[$count][2]++; #1-fold count
			$T_regions{$chr}[$count][3]+=$cov; #total coverage per region (1-fold)
			$totalcov+=$cov; #add up total coverage
			if($cov>$max){$max=$cov;}#set max
			if($cov<$min){$min=$cov;}#set min
			$Startfrom{$chr} = $count; last;
			}
	}#end of for target regions loop
}#end of while pileup file loop
close PILEUP;

my $mean_cov_TT = $totalcov/$Tbase_count;

my $pct_Bases1fold = ($Bcount_1fold/$Tbase_count)*100;
my $pct_Bases5fold = ($Bcount_5fold/$Tbase_count)*100;
my $pct_Bases10fold = ($Bcount_10fold/$Tbase_count)*100;
my $pct_Bases20fold = ($Bcount_20fold/$Tbase_count)*100;
my $pct_Bases30fold = ($Bcount_30fold/$Tbase_count)*100;

my $num_targets_total=0;
my $num_targets=0;
my $num_targets30=0;
my $num_targets20=0;
my $num_targets10=0;
my $num_targets5=0;

#open per target output file
my $outfile3 = $path."/Coverage_perTarget".$ID."_alltargets.txt";
open(OUT3, ">>$outfile3") || die "Cannot open file \"$outfile3\" to write to!\n";
print OUT3 "Patient\tChromosome\tTarget Name\tTarget Start Position\tTarget End Position\tTarget Length\tTotal Coverage for all Target Bases\tMean per Base Coverage\tNo. Target Bases Covered (30-fold)\t% Target Bases Covered (30-fold)\tNo. Target Bases Covered (20-fold)\t% Target Bases Covered (20-fold)\tNo. Target Bases Covered (10-fold)\t% Target Bases Covered (10-fold)\tNo. Target Bases Covered (5-fold)\t% Target Bases Covered (5-fold)\tNo. Target Bases Covered (1-fold)\t% Target Bases Covered (1-fold)\n";

foreach my $ch (sort keys %T_regions){
	for (my $ct=0; $ct < $T_perChr{$ch}; $ct++){
		$num_targets_total++;
		my $Tgt_lnth = ($T_regions{$ch}[$ct][1]+1)-$T_regions{$ch}[$ct][0];
		my $mean_per_base_cov=$T_regions{$ch}[$ct][3]/$Tgt_lnth;
		my $pct_bases_1_fold =($T_regions{$ch}[$ct][2]/$Tgt_lnth)*100;
		my $pct_bases_30_fold =($T_regions{$ch}[$ct][8]/$Tgt_lnth)*100;
		my $pct_bases_20_fold =($T_regions{$ch}[$ct][5]/$Tgt_lnth)*100;
		my $pct_bases_10_fold =($T_regions{$ch}[$ct][7]/$Tgt_lnth)*100;
		my $pct_bases_5_fold =($T_regions{$ch}[$ct][6]/$Tgt_lnth)*100;
		#$T_regions{$ch}[$ct][9][$p1]=$pct_bases_10_fold;
		#print per target info to file - all targets
		print OUT3 "$ID\t$ch\t$T_regions{$ch}[$ct][4]\t$T_regions{$ch}[$ct][0]\t$T_regions{$ch}[$ct][1]\t$Tgt_lnth\t$T_regions{$ch}[$ct][3]\t$mean_per_base_cov\t$T_regions{$ch}[$ct][8]\t$pct_bases_30_fold\t$T_regions{$ch}[$ct][5]\t$pct_bases_20_fold\t$T_regions{$ch}[$ct][7]\t$pct_bases_10_fold\t$T_regions{$ch}[$ct][6]\t$pct_bases_5_fold\t$T_regions{$ch}[$ct][2]\t$pct_bases_1_fold\n";
		
		if($Tgt_lnth == $T_regions{$ch}[$ct][2]){$num_targets++;} #num targets completely covered 1-fold
		if($Tgt_lnth == $T_regions{$ch}[$ct][8]){$num_targets30++;} #num targets completely covered 30-fold
		if($Tgt_lnth == $T_regions{$ch}[$ct][5]){$num_targets20++;} #num targets completely covered 20-fold
		if($Tgt_lnth == $T_regions{$ch}[$ct][7]){$num_targets10++;} #num targets completely covered 10-fold
		if($Tgt_lnth == $T_regions{$ch}[$ct][6]){$num_targets5++;} #num targets completely covered 5-fold
				
	#re-initialise hash coverage values to 0
	$Targetregions{$ch}[$ct][2]=0;
	$Targetregions{$ch}[$ct][3]=0;
	$Targetregions{$ch}[$ct][5]=0;
	$Targetregions{$ch}[$ct][6]=0;
	$Targetregions{$ch}[$ct][7]=0;
	$Targetregions{$ch}[$ct][8]=0;
		}}
close OUT3;		

my $pct_targets_1fold = ($num_targets/$num_targets_total)*100; # % targets completely covered 1-fold
my $pct_targets_30fold = ($num_targets30/$num_targets_total)*100; # % targets completely covered 30-fold
my $pct_targets_20fold = ($num_targets20/$num_targets_total)*100; # % targets completely covered 20-fold
my $pct_targets_10fold = ($num_targets10/$num_targets_total)*100; # % targets completely covered 10-fold
my $pct_targets_5fold = ($num_targets5/$num_targets_total)*100; # % targets completely covered 5-fold

#return coverages
my $results_mean_cov = join('	', $ID, $Tbase_count, $mean_cov_TT, $min, $max);
my $results_bases = join('	', $ID, $Tbase_count, $Bcount_30fold, $pct_Bases30fold, $Bcount_20fold, $pct_Bases20fold, $Bcount_10fold, $pct_Bases10fold, $Bcount_5fold, $pct_Bases5fold, $Bcount_1fold, $pct_Bases1fold);
my $results_targets = join('	', $ID, $num_targets_total, $num_targets30, $pct_targets_30fold, $num_targets20, $pct_targets_20fold, $num_targets10, $pct_targets_10fold, $num_targets5, $pct_targets_5fold, $num_targets, $pct_targets_1fold);
return ($results_bases, $results_targets, $results_mean_cov);
}#end of sub

exit;
