#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $file = shift;

my $seqIO_obj = Bio::SeqIO->new (-file=> $file, -format=>'fasta');
while (my $seqObj = $seqIO_obj->next_seq){
  my $name = $seqObj->id;
  my $out_seqIO_obj = Bio::SeqIO->new (-file=> ">$name.fa", -format=>'fasta');
  $out_seqIO_obj->write_seq($seqObj);
} 
