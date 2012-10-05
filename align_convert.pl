#!/usr/bin/perl -w
use strict;
use Bio::AlignIO;
 
my ($in_file, $in_format, $out_file, $out_format) = @ARGV;

if (!defined $in_file or !defined $in_format or !defined $out_format or !defined $out_file){
  die "align_convert.pl in_file in_format out_file out_format\n";
}


my $in  = Bio::AlignIO->new(-file => $in_file ,
                         -format => $in_format);
my $out = Bio::AlignIO->new(-file => ">$out_file" ,
                         -format => $out_format);
 
while ( my $aln = $in->next_aln ) { 
  $out->write_aln($aln); 
}
