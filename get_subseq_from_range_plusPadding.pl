#!/usr/bin/perl -w
use strict;

use Bio::DB::Fasta;

my $dbfile = shift;
my $range_file = shift;
my $padding = shift;

if (!defined $padding){
  $padding = 0;
}
#my @ranges = @ARGV;


if (!defined $dbfile or !defined $range_file){
  die "Please provide fasta file and chromosome postion with ranges\n
example:
./get_subseq_from_range.pl fastaFile range_File

cat range_File
Chr1:12345..23456 [optional name]
Chr4:56789..69012 [optional name]\n";
}

my $db_obj = Bio::DB::Fasta->new($dbfile);

open IN, $range_file or die "Can't open range_file\n";
while (my $range = <IN>){
  chomp $range;
  my ($ref, $start , $end, $name) = $range =~ /^(\S+):(\d+)\.\.(\d+)\s*(.+)?$/;
  if (!defined $ref or !defined $start or !defined $end){
    die "your range is not in the correct format.
your range(s) : 
$range

correct ranges: 
Chr1:2143432..2144432
Chr3:4343432..5144432
"
  }

  my $seq_obj = $db_obj->get_Seq_by_id($ref);
  if ($seq_obj) {
    my $seq = $seq_obj->seq;
    my $padded_start = $start - $padding;
    my $padded_end = $end + $padding;
    if ( defined $seq ) {
      my $subseq = substr $seq, $padded_start-1, ($padded_end - $padded_start) + 1;
      if ($padding == 0){
        my $header =  "$ref:$start..$end";
        if (defined $name){
          $header = "$name $header";
        } 
        print ">$header\n$subseq\n";
      }else{
        my $header =  "$ref:$padded_start..$padded_end|padding|$padding|sub|$ref:$start..$end";
        if (defined $name){
          $header = "$name $header";
        } 
        print ">$header\n$subseq\n";

      }
    }
  }else{
   warn "error retrieving $range\n";
  }
}
