#!/usr/bin/perl -w
use strict;
use Data::Dumper;
use File::Spec;
## parse_blat.pl blatout_dir > blat_parsed_outsummary_results.txt
my $blatout_dir     = shift;
my %seqs;
my %convert;

while (my $line = <DATA>){
  chomp $line;
  my ($SRR,$SRX,$NAM) = split /\s+/ , $line;
  $convert{$SRR}=$NAM;
}


my %found;
my @blat_files = <$blatout_dir/*blatout>;
  foreach my $blat_file (@blat_files) {
  open OUT, ">$blat_file.out.parsed";
  print OUT "qName\tqLen\ttName\ttLen\tm+mm\tmatches\tmismatches\n";
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
    next unless $line =~ /prime/;
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
   ## all seqs are 35 or 36 bp ling
   next if  ($matches + $mismatches) < $qLen;
   ## throw out if there are too many MM 
   ## maybe this it too many
   next if   $mismatches  > 3; 
   my ($library) = $qName =~ /(\w+)\.\d+/;
   $library =~ s/\s//g;
   my $NAM = $convert{$library};
   ##566|ZM_Tourist_22|methyl|chr1:82026867..82026987|complete|chr1:82026867..82026987|five_prime|chr1:82026842..82026891 
   my ($id,$te,$methyl,$methyl_range,$complete,$te_range,$ends,$pos) = split /\|/ , $tName;
   my $te_id = "$id|$te|$te_range";
   #$found{q}{$te_id}{$library}{$qName}{matches}=$matches;
   #$found{q}{$te_id}{$library}{$qName}{mismatches}=$mismatches;
   #$found{q}{$te_id}{$library}{$qName}{len}=$qLen;
   $found{t}{$te_id}{$NAM}{$ends}++;
   
   print OUT "$qName\t$qLen\t$tName\t$tLen\t",$matches+$mismatches,"\t$matches\t$mismatches\n";
 }

}
print "library\tte\tconfirmedEnds\tfive_count\tthree_count\n";
  foreach
    my $te_id ( sort { ( split( /\|/, $a ) )[0] <=> ( split( /\|/, $b ) )[0] }
    keys %{ $found{t}} )
  {
    foreach my $library ( sort keys %{ $found{t}{$te_id} } ) {
      my $five_count =
        exists $found{t}{$te_id}{$library}{five_prime}? $found{t}{$te_id}{$library}{five_prime}: 0;
      my $three_count =
        exists $found{t}{$te_id}{$library}{three_prime}? $found{t}{$te_id}{$library}{three_prime}: 0;
      my $both_ends =
        ( ( $five_count > 0 ) and ( $three_count > 0 ) )? "two-end": "one-end";
      print "$library\t$te_id\t$both_ends\t$five_count\t$three_count\n";
    }
  }


__DATA__
SRR026733	SRX010812	B73
SRR027083	SRX010812	B73
SRR026735	SRX010813	B97
SRR027081	SRX010813	B97
SRR027082	SRX010813	B97
SRR026737	SRX010814	CML103
SRR027077	SRX010814	CML103
SRR027078	SRX010814	CML103
SRR026743	SRX010815	CML228
SRR027075	SRX010815	CML228
SRR027076	SRX010815	CML228
SRR026745	SRX010816	CML247
SRR027073	SRX010816	CML247
SRR027074	SRX010816	CML247
SRR026747	SRX010817	CML277
SRR027072	SRX010817	CML277
SRR026749	SRX010818	CML322
SRR027071	SRX010818	CML322
SRR026751	SRX010819	CML333
SRR027070	SRX010819	CML333
SRR026753	SRX010820	CML52
SRR026758	SRX010821	CML69
SRR027069	SRX010821	CML69
SRR026759	SRX010822	HP301
SRR027036	SRX010822	HP301
SRR026760	SRX010823	Il14H
SRR026736	SRX010824	Ki11
SRR027180	SRX010824	Ki11
SRR026734	SRX010825	Ki3
SRR027034	SRX010825	Ki3
SRR027035	SRX010825	Ki3
SRR026738	SRX010826	Ky21
SRR027031	SRX010826	Ky21
SRR026739	SRX010827	M162W
SRR027020	SRX010827	M162W
SRR026740	SRX010828	M37W
SRR026997	SRX010828	M37W
SRR026742	SRX010830	Mo18W
SRR026783	SRX010830	Mo18W
SRR026784	SRX010830	Mo18W
SRR026744	SRX010831	Ms71
SRR026772	SRX010831	Ms71
SRR026746	SRX010832	NC350
SRR026771	SRX010832	NC350
SRR026782	SRX010832	NC350
SRR026748	SRX010833	NC358
SRR026767	SRX010833	NC358
SRR026752	SRX010835	Oh7B
SRR026766	SRX010835	Oh7B
SRR026754	SRX010836	P39
SRR026765	SRX010836	P39
SRR026756	SRX010837	Tx303
SRR026764	SRX010837	Tx303
SRR026757	SRX010838	Tzi8
SRR026763	SRX010838	Tzi8
