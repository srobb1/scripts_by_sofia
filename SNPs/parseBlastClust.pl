#!/usr/bin/perl -w

use strict;
my $file = shift;    ## blastclust.out file
open IN, $file;
my $db = shift;      ##fasta file of flankers
if ( !-e "$db.phr" ) {
  `formatdb -i $db -p F -o T`;
}

`mkdir -p alignment`;
`rm alignment/*`;
my $line_count = 0;
while ( my $line = <IN> ) {
  $line_count++;
  chomp $line;
  my @reads = split /\s+/, $line;
  my %reads;
  my $clust_count = @reads;
  if ( $clust_count >= 5 ) {
    foreach my $read (@reads) {

      #pong_5prime.192.1:47:12210:48909:Y:end.trim
      my ($TE) = $read =~ /(pong_\dprime)/;
      $reads{$TE}++;
      `fastacmd -d $db -s $read >> alignment/$line_count.$file.fa`;
    }
    my $te_count;
    foreach my $te ( keys %reads ) {
      my $count = $reads{$te};
      $te_count .= "$te:$count\t";
    }
    $te_count =~ s/\t$//;
    print $clust_count , "\t$te_count\n";
`muscle -clw -in alignment/$line_count.$file.fa -out alignment/$line_count.$file.clw 2> muscle.out`;
  }
}
