#!/usr/bin/perl -w
use strict;
use Getopt::Long;
if ( !defined @ARGV ) {
  &getHelp();
}
my $bowtie2          = 0;
GetOptions(
  'p|parallel:i'         => \$parallel,
  'q|qsub_q:s'           => \$qsub_q,
  'h|help' => \&getHelp,
);

sub getHelp {
  print ' 
usage:
./relocaTE.pl [-t TE_fasta_file][-g chromosome_genome_fasta][-d dir_of_fq][-e short_sample_name][-h] 

options:
-t |--te_fasta		file		fasta containing nucleotide sequences of transposable elements with 
';
  exit 1;
}
