#!/usr/bin/perl
# $Id: haystack.pl,v 1.1 2007/05/08 22:21:10 priesth Exp $
#########################################################
# teragrid_wordfinder	                               #
#########################################################

use strict;
use warnings;
use Statistics::Distributions;
use Carp;
use Getopt::Std;
use Cwd;
use File::Copy;
use Math::Complex;
use Time::HiRes qw( gettimeofday tv_interval time );
use vars qw/ $opt_m $opt_d $opt_c $opt_f $opt_r $opt_b $opt_n/;

#########################################################
# Start Variable declarations	                       #
#########################################################

my ($modelfile, $datafile, $correlation_cut, $fold_cut, $pval_cut, $background_cut, $outfile);
my ($outdir);
my $workingdir 		= Cwd::getcwd();
my $time = time();
#########################################################
# End Variable declarations	                         #
#########################################################

#########################################################
# Start Main body of Program	                        #
#########################################################

 getopts(':m:d:c:f:r:b:n:');

 var_check();

open(MODELS,"<", $modelfile ) || die "Cannot open $modelfile!\n$!\nexiting...\n";
my @models=<MODELS>;
close MODELS;

open (OUT,">",$outfile) || die "Cannot open $outfile!\n$!\nexiting...\n";
open (GENE, "<$datafile") || die "Can't open $datafile !!$!\n";
while (<GENE>) {
	chomp($_);
	next unless $_=~m/\d/;	
	my @data;
	if($_=~m/,/){
		@data=split(/\,/,$_);
	}elsif($_=~m/\t/){
		@data=split(/\t/,$_);
	}else{
		#### error handle
	}

	### Get the variables for model re-scaling and/or filtering ###
	my $probeset=shift @data;

	my @sorted = sort { $a <=> $b } @data;
	my $max = $sorted[$#sorted];
	my $min = $sorted[0];
	my $fold=0;
	if ($min > 0) {
		$fold = $max/$min;
	}elsif ($min <= 0) {
		$fold = $max; #prevents div by zero errors
	}
	my $trigger=0;
	my %result;
	if ($max > $background_cut && $fold >= $fold_cut) {
		foreach my $line (@models){
			chomp $line;
			next unless $line=~m/\d/;	
			my @model;
			if($line=~m/,/){
				@model=split(/,/,$line);
			}elsif($line=~m/\t/){
				@model=split(/\t/,$line);
			}else{
			}
			my $model_name=shift @model;
			die "ERROR - The number of MODEL points and DATA points is different!!!\n" if($#model != $#data);
			my @fitted;
			my ($slope, $intercept) = best_line(\@model, \@data);
			for (my $m=0; $m<=$#model; $m++) {
				push @fitted, $slope*$model[$m]+$intercept;
			}
			my $correlation = correlation(\@model, \@data);
			my $r = $correlation;
			my $tstat = $r * sqrt(scalar(@data) - 2)/sqrt(1-($r**2));
			my $tprob=Statistics::Distributions::tprob (scalar(@data) - 2,$tstat); 
	
			if ($correlation > $correlation_cut && $tprob < $pval_cut) {
				$trigger = 1;
				my $dline="";
				$dline.="DATA\t";
				$dline.="$probeset\t";
				$dline.="$model_name\t";
				$dline.="$correlation\t";
				$dline.="$tstat\t";
				$dline.="$tprob\t";
				for (my $n=0; $n<=$#data; $n++) {
					$dline.=$data[$n]."\t";
				}
				$dline.="\n";
				
				my $mline="";
				$mline.="MODEL\t";
				$mline.="$probeset\t";
				$mline.="$model_name\t";
				$mline.="$correlation\t";
				$mline.="$tstat\t";
				$mline.="$tprob\t";
				for (my $o=0; $o<=$#fitted; $o++) {
					$mline.=$fitted[$o]."\t";
				}
				$mline.="\n";
				if ($result{'c'}){
					if ($correlation>$result{'c'}){
						$result{'m'}=$mline;
						$result{'d'}=$dline;
						$result{'c'}=$correlation;
					}else{
					}
				}else{
					$result{'m'}=$mline;
					$result{'d'}=$dline;
					$result{'c'}=$correlation;
				}
			}else {}
		}
	}
	
	if($result{'c'}){
		print OUT $result{'d'};
#		print OUT $result{'m'};
	}else {}
}
close OUT;
close GENE;




### Subroutines ###

sub correlation {
	my ($array1ref,$array2ref) = @_;
	my ($sum1,$sum2);
	my ($sum1_squared,$sum2_squared);
	foreach (@$array1ref) {$sum1 += $_; $sum1_squared += $_ **2}
	foreach (@$array2ref) {$sum2 += $_; $sum2_squared += $_ **2}
	if((((@$array1ref * $sum1_squared) - ($sum1 ** 2))!=0)&&(((@$array1ref*$sum2_squared)-($sum2**2))!=0)){
	 return (@$array1ref **2)*covariance($array1ref,$array2ref) /sqrt(((@$array1ref * $sum1_squared) - ($sum1 ** 2))*((@$array1ref*$sum2_squared)-($sum2**2)));
	}else{
	 return 0;
	}
}

sub covariance {
	my ($array1ref,$array2ref) = @_;
	my ($i,$result);
	for ($i=0; $i<@$array1ref; $i++) {
	    $result += $array1ref->[$i]*$array2ref->[$i];
	}
	$result /= @$array1ref;
	$result -= mean($array1ref)*mean($array2ref);
}

sub mean {
	my ($arrayref) = @_;
	my $result;
	foreach (@$arrayref) { $result += $_ }
	return $result/@$arrayref;
}

sub best_line {
	my ($array1ref, $array2ref) = @_;
	my ($i, $product, $sum1, $sum2, $sum1_squares, $a, $b);
	for ($i = 0; $i < @$array1ref; $i++) {
	    $product += $array1ref->[$i] * $array2ref->[$i];
	    $sum1 += $array1ref->[$i];
	    $sum1_squares += $array1ref->[$i] ** 2;
	    $sum2 += $array2ref->[$i];
	}
	my $incidental=((@$array1ref * $sum1_squares) - ($sum1 ** 2));
	if($incidental==0){
	$b=0;
	}else{
	$b = ((@$array1ref * $product) - ($sum1 * $sum2)) / $incidental;
	}
	$a = ($sum2 - $b * $sum1) / @$array1ref;
	return ($b, $a);
}
# usage: perl profile_matcher.pl model_file data_file correlation_cutoff fold_cut r-to-t_pvalue_cut background_cut
sub var_check {
  if ($opt_m) {
   $modelfile = $opt_m;
  } else {
   var_error();
  }
  if ($opt_d) {
   $datafile = $opt_d;
  } else {
   var_error();
  }
  if ($opt_c) {
   $correlation_cut = $opt_c;
  } else {
   $correlation_cut = ".8";
  }
  if ($opt_f) {
   $fold_cut = $opt_f;
  } else {
   $fold_cut = "2";
  }
  if ($opt_r) {
   $pval_cut = $opt_r;
  } else {
   $pval_cut = ".05";
  }
  if ($opt_b) {
   $background_cut = $opt_b;
  } else {
   $background_cut = "100";
  }
  if ($opt_n) {
   $outfile = $opt_n;
  } else {
   var_error();
  }
}

#########################################################
# Start of Varriable error Subroutine "var_error"	   #
#########################################################

sub var_error {

  print "\n\n";
  print "  haystack.pl identifies the best model for gene expression data, gene by gene.\n\n";
  print "  Usage:\n";
  print "  haystack.pl -m <model file> -d <data file> -c <correlation cutoff value> -f <fold cut value> -r <r-to-t pvalue cutoff> -b <background cutoff value>\n";
  print "\n\n\n";
  print "  -m   The model datafile.\n";
  print "\n";
  print "  -d   The data file.\n";
  print "\n";
  print "  -c   The correlation cutoff value.\n";
  print "\n";
  print "  -f   The fold cutoff value.\n";
  print "\n";
  print "  -r   The r-to-t p-value cutoff value.\n";
  print "\n";
  print "  -b   The background cutoff value.\n";
  print "\n";
  print "  -n   The name of the outputfile.\n";
  print "\n\n\n";
  exit 0;

}

#########################################################
# End of Variable error Subroutine "var_error"	      #
#########################################################