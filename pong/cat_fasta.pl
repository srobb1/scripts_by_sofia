#!/usr/bin/perl -w
use strict;

my @files = @ARGV;
#~/rice/pong/09112012_A123_0/pong/te_containing_fq/*fasta > A123_0.te_containing_fq.fasta
my $count;
foreach my $file (@files){
  next if $file !~ /fasta$/;
  $count++;
  open IN, $file;
  my $seqCount;
  while (my $line = <IN>){
    $seqCount++;
    chomp $line;
   $line =~ s/^>(.+)/>$count.$seqCount $1/;
   print $line , "\n";
  }
}

