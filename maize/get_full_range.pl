#!/usr/bin/perl -w


#10      MTEC    repeat  1883171 1883325 1350    +       .       class=II;subclass=1;order=TIR;superfamily=[DTM] Mutator;family=Zm00118;type=Type II Transposons/TIR;name=DTM_Zm00118_consensus       88.46031746     85.73389623     81.76228403     10      11      49      1       31.72903226

use strict;
my $file = shift;
my $gff = shift;

my %TEs;
open GFF  , $gff or die "Can't open $gff\n";
while (my $line = <GFF>){
  chomp $line;
  next if $line =~/^#/;
  my @line = split /\t/ , $line;
  my ($ref , $start, $end, $nineth) = ($line[0] , $line[3], $line[4], $line[8]);
  my ($te) = $nineth =~ /name=(\S+)/;
  if (!defined $start or !defined $end){
    print "Weird $line\n";
  }
  $TEs{$te}{$ref}{$start}="$start..$end";
  $TEs{$te}{$ref}{$end}="$start..$end";
}
open IN , $file or die "Can't open $file\n";
my %methyl;
while (my $line = <IN>){
  chomp $line;
  next if $line =~/seqname/;
  my @line = split /\t/ , $line;
  my ($ref , $start, $end, $ninth) = ($line[0] , $line[3], $line[4], $line[8]);
  if (!defined $start or !defined $end){
    print "Weird $line\n";
    die;
  }
  my ($te) = $ninth =~ /name=(\S+)/;
  if (exists $TEs{$te}{$ref}{$start}){
    my $range = $TEs{$te}{$ref}{$start};
    push @{$methyl{$te}{"chr$ref:$range"}} , "chr$ref:$start..$end";
  }elsif(exists $TEs{$te}{$ref}{$end}){
    my $range = $TEs{$te}{$ref}{$end};
    push @{$methyl{$te}{"chr$ref:$range"}} , "chr$ref:$start..$end";
  }else{
    print "Not found $te, $ref, $start, $end\n";
  }
}
my $count = 0;
foreach my $te (keys %methyl){
  foreach my $complete (keys %{$methyl{$te}}){
    my $methyl = join '|', @{$methyl{$te}{$complete}};
    $count++;
    print "$count|$te|methyl|$methyl|complete|$complete\n";
  }
}
