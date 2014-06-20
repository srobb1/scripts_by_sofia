#!/usr/bin/perl -w
use strict;

my @files = @ARGV;
#Chr10	HEG4_2	transposable_element_insertion_site	21716391	21716819	.	.	.	ID=mping.te_insertion_site.Chr10.21716391..21716819;TE_Name=mping;Note=Shared, in ref and reads;left_flanking_read_count=1;right_flanking_read_count=0
my %features;
foreach my $file (@files){
  open GFFIN, $file or die "Can't open $file\n";
  while (my $line =<GFFIN>){
    chomp $line;
    next if $line =~ /^#/;
    ## should only be Shared and reference-only insertions
    next if $line =~ /Non-ref/;
    my ($ref,$source,$feature,$start,$end,$score,$strand,$phase,$nine) = split /\t/ , $line;
    my ($TE) = $nine =~ /ID=([A-Za-z0-9_-]+)/; 
    my $note = $nine=~/Reference-only/ ? 'Reference-only' : 'Shared';
    if (!exists $features{"$TE.te_insertion_site.$ref.$start..$end"}{Shared} and $note eq 'Reference-only'){
      $features{"$TE.te_insertion_site.$ref.$start..$end"}{$note}{"$ref\tNB\ttransposable_element_insertion_site\t$start\t$end\t$score\t$strand\t$phase\tID=$TE.te_insertion_site.$ref.$start..$end;Note=$note"}++;
   }elsif($note eq 'Shared' and exists $features{"$TE.te_insertion_site.$ref.$start..$end"}{'Reference-only'}){
      delete  $features{"$TE.te_insertion_site.$ref.$start..$end"}{'Refenence-only'};
      $features{"$TE.te_insertion_site.$ref.$start..$end"}{$note}{"$ref\tNB\ttransposable_element_insertion_site\t$start\t$end\t$score\t$strand\t$phase\tID=$TE.te_insertion_site.$ref.$start..$end;Note=$note"}++;
   }else {
     $features{"$TE.te_insertion_site.$ref.$start..$end"}{$note}{"$ref\tNB\ttransposable_element_insertion_site\t$start\t$end\t$score\t$strand\t$phase\tID=$TE.te_insertion_site.$ref.$start..$end;Note=$note"}++;
   }
    print $line,"\n" if $line =~ /Shared/;
  }
}

foreach my $feat (keys %features){
  foreach my $type (keys %{$features{$feat}}){
    foreach my $line (keys %{$features{$feat}{$type}}){
      print "$line\n";
    }
  }
}
