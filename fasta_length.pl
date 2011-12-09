#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $file = shift;

my $seqIO_obj = Bio::SeqIO -> new (-file => $file , -format => 'fasta');

my $total_len;
while (my $seq_obj = $seqIO_obj->next_seq){
	my $id = $seq_obj->id;
	my $len = $seq_obj->length;
	$total_len += $len;
	print "$id\t$len\n";

}
print "totalLen\t$total_len\n";
