#!/usr/bin/perl -w
use strict;

my $file = shift;

open IN, $file or die "cant open $file\n";
my $i;
while (my $line = <IN>){
  next unless $line =~ /^primers/;
  chomp $line;
  $i++;
  #primers	LOC_Os01g66970(419)	Chr1:38895496..38895527,38895887..38895914
  my ($type, $name, $range) =split /\s+/ , $line;
  my ($ref,$s1,$e1,$s2,$e2) = $range =~ /(\w+):(\d+)\.\.(\d+),(\d+)\.\.(\d+)/;
  my ($product_size) = $name=~/\((\d+)\)/;
  my $ID="Primer.$i.$ref.$s1.$e2";
  my $C1_ID="Primer.$i.$ref.$s1.$e1";
  my $C2_ID="Primer.$i.$ref.$s2.$e2";
  print join("\t",$ref,'.','match',$s1,$e2,'.','.','.',"ID=$ID;product_size=$product_size"),"\n";
  print join("\t",$ref,'.','match_part',$s1,$e1,'.','.','.',"ID=$C1_ID;Parent=$ID"),"\n";
  print join("\t",$ref,'.','match_part',$s2,$e2,'.','.','.',"ID=$C2_ID;Parent=$ID"),"\n";
}
