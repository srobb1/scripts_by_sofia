#!/usr/bin/perl -w
use strict;
use Data::Dumper;
my $repeatMakerOut = shift; 
my $vcf_tab = shift;
open INRM , "$repeatMakerOut" or die "Can't open RepeatMasker file $repeatMakerOut, $!\n";
open VCFTAB , "$vcf_tab" or die "Can't open RepeatMasker file $vcf_tab, $!\n";

my %repeats;
<INRM>; #throw out header
<INRM>; #throw out header
<INRM>; #throw out header
while (my $line = <INRM>){
  chomp $line;
  my @line = split /\s+/ , $line;
  my ($chr , $start, $end) = ($line[5] , $line[6] , $line[7]);
  $repeats{$chr}{$start}{start} = $start;
  $repeats{$chr}{$start}{end} = $end;
}
#print Dumper \%repeats;
#my $header = <VCFTAB>; #throw out header
#print $header;
my $count = 0;
while (my $line = <VCFTAB>){
  if ($line =~ /^#/){
   print $line ;
   next;
  } 
  chomp $line;
  my @line = split /\t/ , $line;
  my ($chr, $pos) = ($line[0],$line[1]);
  my $keep = 1;
  ##is SNP in a repeat range?
  foreach my $repeat_start (keys %{$repeats{$chr}}){
    #print "$repeats{$chr}{$repeat_start}{start} -- $repeats{$chr}{$repeat_start}{end} \n";
    if ( $pos >= $repeats{$chr}{$repeat_start}{start} and $pos <= $repeats{$chr}{$repeat_start}{end}){
      $keep = 0;
      $count++;
      last;
    }
  }
  print "$line\n" if $keep;
}
open OUT, ">$vcf_tab.repeat_count_thrown_out";
print OUT "threw out $count SNPs in reapeat regions\n";
