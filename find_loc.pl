#!/usr/bin/perl -w 
use strict;

my $file = shift;

open IN, $file;

`cat  0.1MM/07182012_HEG4_2_strand/*/results/Chr*all* > HEG4_all.sites.txt`;
`cat  0.1MM/07188012_EG4_strand/*/results/Chr*all* > EG4_all.sites.txt`;


##Gaijin  Chr10.19836983  1       EG4_1
##Gaijin  Chr10.5025393   2       EG4_1,HEG4_2
while (my $line = <IN> ){
  chomp $line;
  next unless $line =~ /Chr/;
  my ($TE, $pos, $count ,$strain) = split /\t/ , $line;
  next unless $count == 1;
    if ($strain =~ /^HEG4/){
       my $output = `grep $pos EG4_all.sites.txt`;
       if (length $output < 2){
          print "$line\n";#:\n$output\n\n"; 
       }
    }else{
       my $output = `grep $pos HEG4_all.sites.txt`;
       if (length $output < 2){
          print "$line\n";#:\n$output\n\n"; 
       }
       #print "$line:\n$output\n\n"; 
    }
}
