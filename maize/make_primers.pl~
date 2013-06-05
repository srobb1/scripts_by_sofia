#!/usr/bin/perl -w
use strict;
use Bio::SeqIO;

my $fasta = shift;
my $cwd =  `pwd`;
$cwd =~ s/\n//;
if (!-d "$cwd/primers"){
  `mkdir -p $cwd/primers`;
}

my $in = "$cwd/primers/$fasta.primer3in";
open P3IN, ">$in";

my $seqIO_obj = Bio::SeqIO -> new ('-format' => 'fasta' , '-file'=>$fasta);
while (my $seq_obj = $seqIO_obj->next_seq){
  my $name = $seq_obj -> id;
  # >1|ZM_CACTA_72|padding|400|te|chr8:118402924..118403878|pcr|chr8:118402524..118404278
  my ($id,$te,$padding,$loc) = $name =~ /^(\d+)\|(.+)\|padding\|(\d+)\|te\|(.+)/;
  my $seq = $seq_obj -> seq;
  my $len = length $seq;
  my ($start , $end) = ($padding+1 , ($len - $padding));
  print P3IN "SEQUENCE_ID=$id|$te|$loc|1
SEQUENCE=$seq
SEQUENCE_TARGET=$start,$end
PRIMER_OPT_SIZE=22
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_PRODUCT_SIZE_RANGE=200-500 600-$len
PRIMER_GC_CLAMP=1
=\n";
  ##if TE is more than 0.5kb, make internal primers
  if ( ($len - ($padding*2)) > 500){
    ($start , $end) = ($padding/2 , ($padding + ($padding/2)));

    print P3IN "SEQUENCE_ID=$id|$te|$loc|2
SEQUENCE=$seq
SEQUENCE_TARGET=$start,$end
PRIMER_OPT_SIZE=22
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_PRODUCT_SIZE_RANGE=200-500 600-$len
PRIMER_GC_CLAMP=1
=\n";
    ($start , $end) = ($len - ($padding + ($padding/2)) , ($len - ($padding/2)));
    print P3IN "SEQUENCE_ID=$id|$te|$loc|3
SEQUENCE=$seq
SEQUENCE_TARGET=$start,$end
PRIMER_OPT_SIZE=22
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_PRODUCT_SIZE_RANGE=200-500 600-$len
PRIMER_GC_CLAMP=1
=\n";
  }
}

