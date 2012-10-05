#!/usr/bin/perl -w

## this retrieves the sequence of an "empty_site"  TE sequence based on start and stop in excel file

use strict;
use Bio::DB::Fasta;
use Bio::Seq;

my $file  = shift;
my $fasta = shift; ##genome fasta

if ( !defined $file or !defined $fasta ) {
  die "./script file_with_range_info_to_be_collected fasta_file_of_genome";
}

my %ref;

open IN, "$file";
my @path = split /\// , $file;
my $name = pop @path;
my ($pre,$ext) = $name =~ /(\S+)\.(\S+)$/;
my $path = join '/' , @path;
open FAOUT, ">$path/emptySites.$pre.fa";

while ( my $line = <IN> ) {
  next if $line =~ /Not found/;
  chomp $line;
  ##ZM_CACTA_72|methyl|chr8:118403324..118403478|complete|chr8:118403324..118403478
  my ($name,$methyl,$ref,$start,$end) = $line =~/^(\S+)\|(methyl.+)\|complete\|(\S+):(\d+)\.\.(\d+)$/;
  ##create an empty site
  ## xxxxxxxxxxxTETETETETETETETTETETExxxxxxxxxxxxxxxxx
  ## xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  my $subseq_before      = getSeq_fastacmd( $ref, $start-20,$start-1 );
  my $subseq_after      = getSeq_fastacmd( $ref, $end+1, $end+20 );
  print FAOUT ">$name|$methyl|complete|$ref:$start..$end|empty_site\n$subseq_before"."$subseq_after\n";
}

sub getSeq_fastacmd {
  my $ref   = shift;
  my $start = shift;
  my $end   = shift;
  my $seq;
  if (!exists $ref{$ref}){ 
    my $record = `fastacmd -d $fasta -s $ref`;
    my ($header, @seq) = split /\n/ , $record; 
    $seq = join '', @seq;
    $ref{$ref}=$seq;
  }else{
    $seq = $ref{$ref};
  
  }
  if ($seq){
    my $length = length($seq);
    $end = $length < $end ? $length : $end;
    if ( defined $seq ) {
      ##correct for 0 1 base notation
      $start--;
      $end--;
      my $subseq = substr ($seq,$start, ($end-$start+1));
      return $subseq;
    }
    else {
      return "subseq retrival issue";
    }
  }else {
    return "seqObj retrieval issue";
  }

}
