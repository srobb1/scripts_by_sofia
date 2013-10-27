#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;
## takes many single OG fasta files
my @files = @ARGV;

my %concat;
open LENS, ">concat_lens.txt";
open SEQS, ">concat_seqs.fa";
open OGS, ">concat_OGs.txt";
foreach my $file (@files){
  my @file = split /\// , $file;
  my $OG = pop @file;
  my $seqIO_obj = Bio::SeqIO -> new ( -file=> $file , -format=>'fasta');
  while (my $seq_obj = $seqIO_obj -> next_seq){
    my $id = $seq_obj -> id;
    my $seq = $seq_obj -> seq;
    my $len = $seq_obj -> length;
    my ($spe) = $id =~ /(\S+)_prot\|/;
    if ( $id =~ /(\S+)_prot\|/){
      $spe = $1;
    }else {
      $spe = $id;
    }
    $concat{$spe}{$OG}{seq}=$seq;
    $concat{$spe}{$OG}{len}=$len;
    $concat{$spe}{$OG}{count}++;
  }
}
  foreach my $spe ( sort keys %concat){
    my $lengths;
    my $seqs;
    my $OGs;
    foreach my $OG ( sort keys %{$concat{$spe}}){
      my ($OG_id) = $OG =~ /^(\w+)\./;
      my $count = $concat{$spe}{$OG}{count};
      my $len = $concat{$spe}{$OG}{len};
      my $seq = $concat{$spe}{$OG}{seq};
      if ($count == 1){
        $lengths  .= "($len)";
        $seqs  .= "$seq";
        $OGs .= "($OG_id)";
      }else{
        delete  $concat{$spe}{$OG};
      } 
    }
    print LENS ">$spe\n$lengths\n";
    print SEQS ">$spe\n$seqs\n";
    print OGS ">$spe\n$OGs\n";
  }




