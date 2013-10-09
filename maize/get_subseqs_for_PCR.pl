#!/usr/bin/perl -w

## this retrieves the compete TE sequence based on start and stop in excel file

use strict;
use Bio::DB::Fasta;
use Bio::Seq;

my $file  = $ARGV[0];
my $fasta = $ARGV[1]; ##genome fasta
my $padding = defined $ARGV[2] ? $ARGV[2]  : 400;
print $padding ,"\n";
if ( !defined $file or !defined $fasta ) {
  die "./script file_with_range_info_to_be_collected fasta_file_of_genome\n";
}
my %ref;
open IN, "$file";
my @path = split /\// , $file;
my $name = pop @path;
my ($pre,$ext) = $name =~ /(\S+)\.(\S+)$/;
my $path = join '/' , @path;
open FAOUT, ">for_PCR_TEplus$padding"."_seq.$pre.fa";
while ( my $line = <IN> ) {
  next if $line =~ /Not found/;
  chomp $line;
  #ZM_CACTA_72|methyl|chr8:118403324..118403478|complete|chr8:118403324..118403478
  #my ($name,$methyl,$ref,$te_start,$te_end) = $line =~/^(\S+)\|(methyl.+)\|complete\|(\S+):(\d+)\.\.(\d+)$/;
  #>10213_ZM_mhAT_288 TSD=........ chr2:20658342..20659773,len=1432,score=1.9,strand=+,methylCG=95.30669146,methylCHG=71.24142661,methylCHH=0.239051444
  my ($name,$ref,$te_start,$te_end) = $line =~/^>(\S+).+(chr\d+):(\d+)\.\.(\d+)/;
  my $start = $te_start - $padding;
  my $end = $te_end + $padding;
  my $subseq      = getSeq_fastacmd( $ref, $start,$end );
  print FAOUT ">$name|padding|$padding|te|$ref:$te_start..$te_end|pcr|$ref:$start..$end\n$subseq\n";
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
