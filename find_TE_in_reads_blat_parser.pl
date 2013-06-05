#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Spec;
## parse_blat.pl blatout_dir > blat_parsed_outsummary_results.txt

my $blat_file = shift;
my %TE;

##blat parser
open INBLAT, $blat_file, or die "Please provide a blat output file\n";

<INBLAT>;    #get rid of blat header lines
<INBLAT>;
<INBLAT>;
<INBLAT>;
<INBLAT>;
while ( my $line = <INBLAT> ) {

  chomp $line;
  next if $line =~ /Searched|Loaded/;
  next if $line =~/^\s*$/;
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
  next if ( $matches + $mismatches ) < $qLen;
  ## throw out if there are too many MM
  ## maybe this it too many
  next if $mismatches > 3;
  ## throw out if the gap is too big
  next if $qBaseInsert > 3;
  next if $tBaseInsert > 3;

  if ( !exists $TE{$tName} ) {
    ${ $TE{$tName}{coverage} }[0] = $tLen;
    for ( my $i = 1 ; $i < $tLen + 1 ; $i++ ) {
      ${ $TE{$tName}{coverage} }[$i] = 0;
    }
  }
  ##if it passes above tests, add 1 for read that algins to TE
  $TE{$tName}{count}++;
  for ( my $i = $tStart ; $i < $tStart + 1 ; $i++ ) {
    ${ $TE{$tName}{coverage} }[$i] = 1;
  }
}
print "TE\tAlignedNT\tTE_len\t%coverage\tReadCount\n";
foreach my $TE ( sort keys %TE ) {
  my $len     = shift @{ $TE{$TE}{coverage} };
  my $matches = 0;
  foreach my $pos ( @{ $TE{$TE}{coverage} } ) {
    $matches += $pos;
  }
  my $count = $TE{$TE}{count};
  my $coverage = $matches / $len;
  print "$TE\t$matches\t$len\t$coverage\t$count\n";
}
