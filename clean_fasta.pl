#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
my $fasta = shift;
die "Please provide a fasta file to be cleaned up\n" if !defined $fasta;

my $seqIO_obj = Bio::SeqIO -> new (-file => $fasta , -format => 'fasta');
my $seqIO_out_obj = Bio::SeqIO -> new (-file => ">$fasta.clean" , -format => 'fasta');
while (my $seq_obj = $seqIO_obj -> next_seq){
  #my $id = $seq_obj -> id;
  #my $desc = $seq_obj -> desc;
  #my $seq = $seq_obj -> seq;
  #print ">$id $desc\n$seq\n";
  $seqIO_out_obj->write_seq($seq_obj);
}
