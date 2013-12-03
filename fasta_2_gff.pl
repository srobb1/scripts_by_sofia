#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $file = shift;

my $seqIO_obj = Bio::SeqIO -> new (-file => $file , -format => 'fasta');

my $total_len;
print "##gff version 3\n";
while (my $seq_obj = $seqIO_obj->next_seq){
	my $id = $seq_obj->id;
        my $desc = $seq_obj->desc;
        my $note = defined $desc ? "Note=$desc" : "";
	my $len = $seq_obj->length;
	$total_len += $len;
  
	print "$id\t.\tchromosome\t1\t$len\t.\t.\t.\tID=$id;Name=$id;Note=$note;\n";
}
