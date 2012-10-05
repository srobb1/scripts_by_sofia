#!/usr/bin/perl -w 
use strict;

my $bam = shift;
my $ref_fasta = shift;


my %ref;
open FASTA,$ref_fasta;
my $count =0;
my $seqCount=0;
my $seq;
while (my $line = <FASTA>){
  if ($line =~ /^>/){
    $seqCount++;
    last if $seqCount > 1;
    next;
  }
  chomp $line;
  $seq .= $line;
}
my @seq = split '' , $seq;
foreach my $nt (@seq){
  $count++;
  $ref{$count}=$nt;
}



my @sam_lines = `samtools view $bam`;
my %ping;
foreach my $sam_line (@sam_lines){
  #1:44:14138:11631:Y      1024    ping    1       255     76M     *       0       0       GGCCAGTCACAATGGGGGTTTCACTGGTGTGTCATGCACATTTAATAGGGGTAAGACTGAATAAAAAATGATTATT    HHHHHHBHHHGBEEGHHHHHHGHHHDE8D>DBB@DFHHHHHHEHGHHGG@8>C?C@CFFE>FFFFACCC@EBABEB    XA:i:1  MD:Z:15A60     RG:Z:HEG4       NM:i:1
  next if $sam_line =~ /^@/;
  my ($read,$flag,$ref,$start,$score,$cigar,$one,$two,$three,$seq,$qual) = split /\t/ , $sam_line;
  next if $ref eq '*';
  my @seq = split '' , $seq;
  my @qual = split '' , $qual;
  for (my $i = 0 ; $i < @seq ; $i++){
    my $nt = $seq[$i];
    my $nt_qual = $qual[$i];
    my @q_scores = ( '!' , "\"" , '#' , '$');# , '%' , '&' , "\'" );
    my $bad_base = 0;
    foreach my $q_score ( @q_scores ){
      $bad_base = 1 if $nt_qual eq $q_score ;
      last if $bad_base;
    }
    next if $bad_base;
    $ping{$start}{$nt}++;
    $start++;
  }
}
my @nts = qw (A T G C N);
foreach my $start (sort {$a <=> $b} keys %ping){
  print "$start($ref{$start})\t";
  my $depth;
  foreach my $nt (@nts){
    my $count = exists $ping{$start}{$nt} ? $ping{$start}{$nt} : 0;
    print "$nt:$count\t";
    $depth+=$count;
  }
  print "$depth\n";
}
