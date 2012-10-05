#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

## each sequence is 20bp before the TE and 20bp after TE, for a total of 40 bp (20/20).
## --------------------TETETETETETETETETETETETET--------------------
## --------------------|-------------------
##
## each sequence is compared to every other sequence, but only 20bp of it, 10bp before the TE and 10bp after the TE (10/10)
##
##           ----------|----------
##
## only the seqs which are uniqe at the 10/10 bp level are kept. this is to insure uniq mapping of reads. reads are 35bp long
## and are required to map at least 10 of the before and after seq


my $fasta = shift;
#>250|ZM_CACTA_72|empty|chr1:7296691..7296740
#ACAAACGAACAAGCCTAAGGGGCCGTTCGTTTGTTTCGGAATGGAGGCCT

my %seqs;

my $seqIO_obj = Bio::SeqIO -> new (-file=> $fasta , -format=>'fasta');
open DUPS, ">dups.$fasta.txt";
open UNIQ, ">uniq.$fasta";
while (my $seq_obj = $seqIO_obj -> next_seq){
  my $id = $seq_obj->id;
  my $subseq = $seq_obj->subseq(10,30);
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
    print DUPS "$dups\n";
  }
}
