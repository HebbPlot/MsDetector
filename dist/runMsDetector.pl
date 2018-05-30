#! /usr/bin/perl -w

# A script to run MsDetector

use strict;

#
# Please, specify the following parameters
#

my $hmm = ""; # HMM file
my $glm = ""; # GLM file
my $seq = ""; # input sequence file
my $msk = ""; # output masked sequence file
my $rpt = ""; # output MSs file
my $scr = ""; # output scores file
my $msDetector = '/path/MsDetector'; # executable
my $len = 6; # the length of the motif
my $fct = 4; # window size = $len x $fct
my $mtr = 'Id'; # scoring matrix
my $frq = "";   # the nucleotide probabilities file

my $cmd = "$msDetector -seq $seq ";
$cmd = $cmd . "-hmm $hmm ";
$cmd = $cmd . "-glm $glm ";
$cmd = $cmd . "-len $len ";
$cmd = $cmd . "-fct $fct ";
$cmd = $cmd . "-mtr $mtr ";
$cmd = $cmd . "-msk $msk ";
$cmd = $cmd . "-rpt $rpt ";
$cmd = $cmd . "-thr 0.5 ";
$cmd = $cmd . "-scr $scr ";
$cmd = $cmd . "-frq $frq ";

print "$cmd\n";

if ( system($cmd) ) {
	print "Could not execute $cmd\n";
}
