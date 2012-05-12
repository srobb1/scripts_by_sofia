#!/usr/bin/perl -w
use strict;

use Bio::SeqIO;

my $file = shift;
my $start =shift;
my $end = shift;

my $seqIO = Bio::SeqIO -> new (-file => "$file" , -format => 'fasta');
while (my $seq_obj = $seqIO -> next_seq){
 my $id = $seq_obj->id;
 my $subseq = $seq_obj->subseq($start,$end);
 print ">$id:$start..$end\n$subseq\n";
}
