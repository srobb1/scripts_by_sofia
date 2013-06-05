#!/usr/bin/perl -w
use strict;
my $dir = shift;
my $cmd = 'tail -n1 '.$dir.'/*process_raw_reads*sh.e* | grep \'\[\' |  sort | uniq -c';
my @out = `$cmd`;
foreach my $line (@out){
  chomp $line;
  if ($line =~ /\[bam_sort_core\] merging from/ or
    $line =~ /\[samopen] SAM header is present:/){
    next;
  }else{
    die "was not expecting: $line\n";
  }
}
