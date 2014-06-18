#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $fasta = shift;  ## 600nt - tsd - mping - tsd - 600nt
my $tsd_len = 3;
my $padding = 600 - $tsd_len;
my $TE_len = 430;

my $seqIO_obj = Bio::SeqIO->new(-file=>$fasta, -format=>'fasta');
while (my $seq_obj = $seqIO_obj->next_seq){
  my $id = $seq_obj->id;
  my $desc = $seq_obj->desc;
  my $seq = $seq_obj->seq;
  #print "($id,$desc,$seq)\n";
  #print substr($seq,0,597),"\n";
  #print substr($seq,597,3),"\n";
  #print substr($seq,600,430),"\n";
  #print substr($seq,1030,3),"\n";
  #print length(substr($seq,1033)),"\n";
  #my ($F1,$TSD,$F2) = $seq =~ /(.{10})(.{10}).{10}.{10}(.{10})/;
  my ($F1,$TSD,$F2) = $seq =~ /(.{$padding})(.{$tsd_len}).{$TE_len}.{$tsd_len}(.{$padding})/;
  #print "$id\t$F1\t$TSD\t$F2\n";
  print ">empty.$id\n",lc($F1),uc($TSD),lc($F2),"\n";
  
}
