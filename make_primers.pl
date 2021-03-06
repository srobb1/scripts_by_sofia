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
  #10811_ZM_mhAT_72|padding|400|te|chr2:43338976..43339669|pcr|chr2:43338576..43340069
  # >1|ZM_CACTA_72|padding|400|te|chr8:118402924..118403878|pcr|chr8:118402524..118404278
  my ($id,$te,$padding,$loc) = $name =~ /^(\d+)[|_](.+)\|padding\|(\d+)\|te\|(.+)/;
  my $seq = $seq_obj -> seq;
  my $len = length $seq;
  my ($start , $te_len) = ($padding+1 , ($len - $padding*2));
  my @ranges;
  for (my $i=$te_len+200 ; $i<$len ; $i+=200){
    my $stop = $i+200;
    my $range = "$i-$stop";
    push @ranges, $range;
  }
  #push @ranges, $len-200 . "-$len";
  my $ranges = join " ", @ranges;
  
  print P3IN "SEQUENCE_ID=$id|$te|$loc
SEQUENCE=$seq
SEQUENCE_TARGET=$start,$te_len
PRIMER_OPT_SIZE=22
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=25
PRIMER_OPT_TM=60
PRIMER_MAX_TM=63
PRIMER_MIN_TM=57
PRIMER_OPT_GC=50
PRIMER_MAX_GC=80
PRIMER_MIN_GC=40
PRIMER_PRODUCT_SIZE_RANGE=$ranges
PRIMER_GC_CLAMP=1
=\n";
}

`primer3_core < $cwd/primers/$fasta.primer3in > $cwd/primers/$fasta.primer3out`; 

