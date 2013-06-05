#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Spec;
## parse_blat.pl blatout_dir > blat_parsed_outsummary_results.txt
my $blatout_dir     = shift;
my %found;


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
   next if $mismatches  > 0 ; #3; 
   $found{$tName}{count}++;   
   $found{$tName}{$qName}{qLen}=$qLen;   
   $found{$tName}{$qName}{tLen}=$tLen;   
   $found{$tName}{$qName}{aLen}=$matches + $mismatches;
   $found{$tName}{$qName}{matches}=$matches;   
   $found{$tName}{$qName}{mm}=$mismatches;   
   
 }
}
#print Dumper \%found;
print  "qName\tqLen\ttName\ttLen\treadCount\tm+mm\tmatches\tmismatches\n";
foreach my $tName (sort keys %found){
foreach my $qName (sort keys %{$found{$tName}}){
   next if $qName eq 'count';
   my $count = $found{$tName}{count};
   next if $count < 3;
   my $tLen = $found{$tName}{$qName}{tLen};
   my $qLen = $found{$tName}{$qName}{qLen};     
   my $aLen = $found{$tName}{$qName}{aLen};
   my $matches = $found{$tName}{$qName}{matches};
   my $mm = $found{$tName}{$qName}{mm};  
   print "$qName\t$qLen\t$tName\t$tLen\t$count\t$aLen\t$matches\t$mm\n";
}}
