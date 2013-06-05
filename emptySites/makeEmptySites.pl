#!/usr/bin/perl -w

## this retrieves the sequence of an "empty_site"  TE sequence based on start and stop and TSD length

use strict;
use Bio::DB::Fasta;
use Bio::Seq;

my $file  = shift;
my $TSD_len = shift;
my $padding = shift;
my $fasta = shift; ##genome fasta
my $db_obj = Bio::DB::Fasta->new($fasta);
if ( !defined $file or !defined $fasta or !defined $TSD_len or !defined $padding) {
  die "./script FILE(ref:start..end) TSD_len padding fasta_file_of_genome";
}

my %ref;

open IN, "$file";

while ( my $line = <IN> ) {
  next if $line =~ /TE/;
  chomp $line;
  ##Chr12:2734541..2734970
  my ($te,$ref,$start,$end) = $line =~/(\S+)\s+(\S+):(\d+)\.\.(\d+)/;
  ##create an empty site
  ## xxxxxxxxxxxTETETETETETETETTETETExxxxxxxxxxxxxxxxx
  ## xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  my $empty_1 = $start-$padding;
  my $empty_2 = $start -1;
  my $empty_3 = $end+$TSD_len+1;
  my $empty_4 = $empty_3+$padding-1;
  my $subseq_before      = getSeq_fastacmd( $ref, $empty_1 , $empty_2 );
  my $subseq_after      = getSeq_fastacmd( $ref, $empty_3 , $empty_4 );
  print ">$ref:$empty_1..$empty_2","[$te($ref:$start..$end)]$empty_3..$empty_4\n$subseq_before"."$subseq_after\n";
}

sub getSeq_fastacmd {
  my $ref   = shift;
  my $start = shift;
  my $end   = shift;
  my $seq_obj= $db_obj->get_Seq_by_id($ref);
  if ($seq_obj){
    my $length = $seq_obj->length;
    $end = $length < $end ? $length : $end;
    if ( defined $seq_obj ) {
      my $subseq  = $seq_obj->subseq($start=>$end);
      return $subseq;
    }
    else {
      return "subseq retrival issue";
    }
  }
}
