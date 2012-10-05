#!/usr/bin/perl -w
use strict;
my $file = shift;
my $before = shift;
my $after = shift;

if (!defined $before){
  $before = 8;
}
if (!defined $after){
  $after = 10;
}
open IN, $file;
while (my $line = <IN>){
  chomp $line;
  my $pos;
  if ($line =~ /^>.+(Chr\d+\.\d+)/){
    $pos = $1;
    chomp (my $seq = <IN>);
    my $mid;
    my $TSD;
    if ($seq =~ /.{$before}(...).{$after}/){
      $mid =  $1;
      if ($mid eq 'TAA'){
        $seq =~ tr/ATGCN/TACGN/;
        $seq = reverse $seq;
        my $subseq = substr $seq, length $mid;
        $seq = $subseq;
        print ">$pos revcomp TSD=$mid\n$seq\n";
      }elsif ($mid eq 'TTA'){
        $TSD = $mid;
        #my $subseq = substr $seq, 0, ((length $seq) -2);
        my $subseq = substr $seq, 0;
        $seq = $subseq;
        print ">$pos TSD=$mid\n$seq\n";
      }else {
        print "something diff $mid $seq\n";
      }
      
    }
  }
}
