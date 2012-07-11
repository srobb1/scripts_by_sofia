#!/usr/bin/perl -w

use strict;
use Bio::DB::Fasta;

my $in_fa       = shift;
my $padding     = shift;
my @coordinates = @ARGV;
if ( !defined $padding or  @coordinates < 1 or !defined $in_fa ) {
  die "Please provide a fastaFile, cooridinate (ex. Chr1.28184849) and a Padding
example:
	get_subseq_from_postition_BATCH.pl Rice.fa 700 Chr2.18713620 Chr3.13232706 Chr10.2854155 Chr1.39505222 Chr1.6172235 Chr10.18001938 Chr3.35945981 Chr7.29388014 Chr3.25586159 Chr3.28019802 Chr11.24682794 Chr1.38092040 Chr10.10427403

";
}
my $dbfile = "$in_fa";
my $db_obj = Bio::DB::Fasta->new($dbfile);

foreach my $coordinate ( sort @coordinates ) {
  my ( $ref, $position ) = $coordinate =~ /(\S+)\.(\d+)/;

  # retrieve a sequence
  my $seq_obj = $db_obj->get_Seq_by_id($ref);
  if ($seq_obj) {
    my $seq = $seq_obj->seq;
    if ( defined $seq ) {
      my $start = $position - $padding - 1;
      my $subseq = substr $seq, $start, ( $padding * 2 ) + 1;
      print ">$coordinate $ref:", $start + 1, "..", $position + $padding + 1,
        "\n$subseq\n";
    }

  }
  else {
    warn("Cannot find $ref\n");
  }
}
