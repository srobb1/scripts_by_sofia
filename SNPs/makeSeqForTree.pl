#!/usr/bin/perl -w
use strict;

my $file =shift;
my $ref = shift;
if (!defined $ref){
  $ref = 'NB';
}
## #CHROM  POS     REF     A119_1  A123_1  EG4_2   HEG4_2
## Chr1    31071   A       G/G     G/G     ./.     G/G

open IN, $file;
my %seqs;
chomp (my $header = <IN>);
my ($head_1, $head_2, $head_3, @strains) =split /\t/, $header;
while (my $line = <IN>){
  next if $line =~/^#/;
  next if $line =~/\.\/\./; ##skip any locations with no data
  chomp $line;
  my ($chr,$pos,$ref_nt,@strain_nt) = split /\t/, $line;
  my @temp;
  my $count=0;
  for (my $i =0 ; $i < @strains ; $i++){
    my ($a1 , $a2) = $strain_nt[$i] =~ /(.)\/(.)/;
    next if $a1 ne $a2;
    $count++;
    $temp[$i] = $a1;    
  }
  if ($count == @strains){
    $seqs{$ref}.=$ref_nt;
    for (my $i =0 ; $i < @strains ; $i++){
      $seqs{$strains[$i]}.= $temp[$i];
    }
  }
  
}
foreach my $strain (sort keys %seqs){
  my $seq = $seqs{$strain};
  my $len = length $seq;
  print ">$strain $len\n$seq\n";
}
