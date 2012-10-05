#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;

my @files = @ARGV;
foreach my $file (@files){
  my ($ref) = $file =~ /\/(.+)\.fasta/;
  my $seqIO = Bio::SeqIO -> new (-file => "$file" , -format => 'fasta');
  while (my $seq_obj = $seqIO -> next_seq){
    my $id = $seq_obj->id;
    my $seq = $seq_obj->seq;
    print ">$ref\t$id\n$seq\n";
  }

}
