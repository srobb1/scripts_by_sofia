#!/usr/bin/perl -w
use strict;

my $tab_file = shift;

open IN, $tab_file;

while (my $line = <IN> ){
  chomp $line;
  if ($line =~ /^#/){
    print "$line\n";
    next; 
  }
  my ($chr, $pos, $ref , @strains) = split /\t/ , $line;
  my $het = 0;
  foreach my $strain (@strains){
    my ($a1 , $a2) = $strain =~ /(.)\/(.)/;
    if ($a1 ne $a2){
      $het = 1;
      last;
    }
  } 
  print "$line\n" if !$het;
}

