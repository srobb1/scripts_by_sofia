#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $file = shift;
#open OUTFASTA , ">new.fa";
my $seqIO_obj = Bio::SeqIO->new (-file=> $file, -format=>'fasta');
my $out_seqIO_obj = Bio::SeqIO->new (-file=> ">new.fa", -format=>'fasta');
while (my $seqObj = $seqIO_obj->next_seq){
  my $len = $seqObj->length;
  my $name = $seqObj->id;
  my $desc = $seqObj->description;
  my $seq = $seqObj->seq;
  next if $len == 0;
  next if $desc =~ /whole genome shotgun sequence/;
  next if $desc =~ /BAC/;
  next if $desc =~ /SEQUENCING IN PROGRESS/;
  next if $desc =~ /Aedes aegypti, clone .+, complete sequence/;
#  print OUTFASTA ">$name $desc\n$seq\n";
  $out_seqIO_obj->write_seq($seqObj) if $len != 0;

} 
