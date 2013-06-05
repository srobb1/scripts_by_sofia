#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Spec;
## parse_blat.pl blatout_dir > blat_parsed_outsummary_results.txt
my $blatout_dir     = shift;
my %found;


  print  "qName\tqLen\ttName\ttLen\tm+mm\tmatches\tmismatches\n";
my @blat_files = <$blatout_dir/*blatout>;
  foreach my $blat_file (@blat_files) {
#  open OUT, ">$blat_file.out.parsed";
#  print OUT "qName\tqLen\ttName\ttLen\tm+mm\tmatches\tmismatches\n";
  my @file_path = split '/', $blat_file;
  my $file_name = pop @file_path;
  ##blat parser
  open INBLAT, $blat_file, or die "Please provide a blat output file\n";

  <INBLAT>;    #get rid of blat header lines
  <INBLAT>;
  <INBLAT>;
  <INBLAT>;
  <INBLAT>;
  while ( my $line = <INBLAT> ) {
    my @line        = split /\t/, $line;
    my $matches     = $line[0];
    my $mismatches  = $line[1];
    my $qBaseInsert = $line[5];
    my $tBaseInsert = $line[7];
    my $strand      = $line[8];
    my $qName       = $line[9];
    my $qLen        = $line[10];
    my $tName       = $line[13];
    my $tLen        = $line[14];
    my $tStart      = $line[15] + 1;
    my $tEnd        = $line[16];
    my $aln_bp      = $matches + $qBaseInsert + $mismatches;
    
   ## throw out if alignment is too small 
   ## blat db is empty sites, query are reads
   next if  ($matches + $mismatches) < $tLen;
   next if $qBaseInsert > 0;
   next if $tBaseInsert > 0;
   ## throw out if there are too many MM 
   ## maybe this it too many
   next if $mismatches  > 3 ;#0 ; #3; 
   $found{$tName}++;   
   print "$qName\t$qLen\t$tName\t$tLen\t",$matches+$mismatches,"\t$matches\t$mismatches\n";
 }
}
print Dumper \%found;
