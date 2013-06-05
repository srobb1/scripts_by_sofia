#!/usr/bin/perl -w
#file:convert_genbank2fasta.pl

use strict;
use Bio::SeqIO;

## ./convert.pl fasta embl infast.fa out.embl

my ($informat,$outformat) = ($ARGV[0],$ARGV[1]);
my ($infile,$outfile) = ($ARGV[2],$ARGV[3]);

my $in = Bio::SeqIO->new(
                -format => $informat,
                -file => $infile,
                );
my $out = Bio::SeqIO->new(
                -format => $outformat,
                -file => ">$outfile"
                );

while ( my $seqObj = $in->next_seq ) {
        $out->write_seq($seqObj);
}

