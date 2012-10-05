#!/usr/bin/perl -w
use strict;

use Bio::DB::Fasta;

my $dbfile = shift;
my @ranges = @ARGV;


if (!defined $dbfile or scalar @ranges eq 0){
  die "Please provide fasta file and chromosome postion with ranges\n
example:
./get_subseq_from_range.pl fastaFile Chr1:12345..23456 Chr4:56789..69012\n";
}

my $db_obj = Bio::DB::Fasta->new($dbfile);
foreach my $range (@ranges){
  my ($ref, $start , $end) = $range =~ /^(\S+):(\d+)\.\.(\d+)$/;
  if (!defined $ref or !defined $start or !defined $end){
    die "your range is not in the correct format.
your range(s) : @ranges
correct ranges: Chr1:2143432..2144432 Chr3:4343432..5144432
"
  }

  my $seq_obj = $db_obj->get_Seq_by_id($ref);
  if ($seq_obj) {
    my $seq = $seq_obj->seq;
    if ( defined $seq ) {
      my $subseq = substr $seq, $start-1, ($end - $start) + 1;
      print ">$ref:$start..$end\n$subseq\n";
    }
  }else{
   warn "error retrieving $range\n";
  }
}
