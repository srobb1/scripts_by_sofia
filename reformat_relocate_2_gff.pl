#!/usr/bin/perl -w
use strict;

#Chr1	A119_2	transposable_element_insertion_site	5648229	5648231	.	-	.	ID=mping.te_insertion_site.Chr1.5648231;Note=Non-reference, not found in reference;left_flanking_read_count=40;right_flanking_read_count=33;left_flanking_seq=ACCAATGAACAATGCCCCTGCATATTCTATCCCCATGTGAATAATTTGATAATTTGATCAATCAGATATCCTAGGATGCATGCCCTGACTTTTCTGCTAA;right_flanking_seq=GTAGGTTGCCACTTGCCCACTACATCTTGGCTTTTTGCCTTCTTCATTATTGCTTTCTCTCTTTTTTGAAAACTTTTCTTCAGTATTGCTTTTGCTATCC;TSD=TTA
#ping    ...     HEG4/EG4        Chr1:4220010..4220012   +

my $file = shift;
open IN , $file or die "Can't open $file $!\n";
while (my $line =<IN>){
  chomp $line;
  my ($te,$tsd,$strains,$loc,$strand) = split /\t/, $line;
  my ($ref,$start,$end) = $loc =~ /(Chr\d+):(\d+)\.\.(\d+)/;
  foreach my $strain (split /\//, $strains){
    print "$ref\t$strain\ttransposable_element_insertion_site\t$start\t$end\t.\t$strand\t.\tID=mping.te_insertion_site.$ref.$end\n";
  }
}
