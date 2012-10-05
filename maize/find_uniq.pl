#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

## each sequence is 25bp outside TE and 25bp inside TE, for a total of 50 bp (25/25).
## -------------------------TETETETETETETETETETETETET
##
## each sequence is compared to every other sequence, but only 20bp of it, 10bp outside of TE and 10bp inside the TE (10/10)
##
##                ----------TETETETETE
##
## only the seqs which are uniqe at the 10/10 bp level are kept. this is to insure uniq mapping of reads. reads are 35bp long
## and are required to map at least 10 to the flanking region and at least 10 bp to the TE of the 25/25 sequence


my $fasta = shift;
#>250|ZM_CACTA_72|five_prime|chr1:7296691..7296740
#ACAAACGAACAAGCCTAAGGGGCCGTTCGTTTGTTTCGGAATGGAGGCCT

my %seqs;

my $seqIO_obj = Bio::SeqIO -> new (-file=> $fasta , -format=>'fasta');
open DUPS, ">dups.$fasta.txt";
open UNIQ, ">uniq.$fasta";
while (my $seq_obj = $seqIO_obj -> next_seq){
  my $id = $seq_obj->id;
  my $subseq = $seq_obj->subseq(15,35);
  $seqs{$subseq}{count}++;
  $seqs{$subseq}{id}{$id}=$seq_obj->seq;
}
#print DUPS "TE\tcount of exact copies of 10/10 region\n";
foreach my $subseq ( keys %seqs){
  my $count = $seqs{$subseq}{count};
  if ($count == 1){
    my ($id) = keys %{$seqs{$subseq}{id}};
    my $seq = $seqs{$subseq}{id}{$id};;
    print UNIQ ">$id\n$seq\n";
  }else{
    my @dups = keys %{$seqs{$subseq}{id}};
    my $dups = join "\n" , @dups;
    #print DUPS "$dups\t$count\n";
    print DUPS "$dups\n";
  }
}
