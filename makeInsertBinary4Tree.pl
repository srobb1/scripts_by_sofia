#!/usr/bin/perl -w
use strict;
#input: file id file id file id ..
#file is output from characterizer, inserts.characterized.txt 
#A123_0  mping   TTA     Chr1.16327646   19.5    0       homozygous
#don't forget to add entries or make a file for NB and the existing TEs in other strains
my %args = @ARGV;
my @ids = values %args;
my %inserts;
foreach my $file (keys %args){
  my $id = $args{$file};
  open IN, $file;
  while (my $line = <IN>){
    chomp $line;
    next if $line =~ /TSD/;
    my ($strain, $te, $tsd, $pos, $flank, $span, $class) = split /\t/ , $line;
    if ($class =~ /homozygous/){
      $inserts{$pos}{$id}=1;
    }else{
      $inserts{$pos}{$id}=0;
    }
    foreach my $other_id(@ids){
      next if $other_id eq $id;
      if (!exists $inserts{$pos}{$other_id}){
         $inserts{$pos}{$other_id}=0;
      }
    }
  }
}
my %seq;
foreach my $pos (sort keys %inserts){
  foreach my $id (keys %{$inserts{$pos}}){
    my $binary = $inserts{$pos}{$id};
    $seq{$id}.=$binary;
  } 
}
foreach my $id (keys %seq){
  my $seq = $seq{$id};
  my $len = length $seq;
  print ">$id $len\n$seq\n";

}
